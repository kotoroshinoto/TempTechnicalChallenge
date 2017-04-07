import sys
import click
import vcf
from exac_interface import make_keystring_for_vcfentry, ExAC, ExAC_VariantData
from vcf.model import _Record as vcfrecord
from typing import List

output_handle = None
exac_vars = list()
exac_defaults = dict()


def add_exac_var(varstr, defval=None):
	exac_vars.append(varstr)
	if defval is not None:
		exac_defaults[varstr] = defval
#by default we want the allele frequency
add_exac_var('variant.allele_freq', 0.0)


def add_exac_var_callback(ctx, param, value):
	print("param: %s\tvalue:%s" % (param.name, value))


class AnnotationRecord:
	def __init__(self, varstring:str):
		self.varkey = varstring
		self.depth = 0
		self.ref_support = 0
		self.ref_support_percent = 0.0
		self.var_support = 0
		self.var_support_percent = 0.0
		self.exac_data = None
		# TODO decide if any additional information available from EXAC is relevant

	def __str__(self):
		cols = [self.varkey.replace('-',"\t"), self.depth,
		        self.ref_support, self.ref_support_percent,
		        self.var_support, self.var_support_percent]
		if self.exac_data is not None:
			cols.append(self.exac_data)
		return "\t".join(str(x) for x in cols)

	@classmethod
	def from_vcf_entry(cls, entry: vcfrecord):
		keystr_list = make_keystring_for_vcfentry(entry)
		records = []
		# freqs = ExAC.get_variant_frequencies(keystr_list)
		for i in range(len(entry.ALT)):
			keystr = keystr_list[i]
			record = cls(keystr)
			# print(keystr, file=output_handle)
			# print(entry.INFO, file=output_handle)
			# print(entry.samples, file=output_handle)
			record.depth = entry.INFO['DP']
			record.ref_support = entry.INFO['RO']
			record.ref_support_percent = float(record.ref_support) / float(record.depth)
			record.var_support = entry.INFO['AO'][i]
			record.var_support_percent = float(record.var_support) / float(record.depth)
			# record.var_exac_frequency = freqs[i]
			records.append(record)
		# TODO poll ExAC for information
		return records


@click.command()
@click.argument('filename', type=click.File('r'))
@click.argument('output', type=click.File('w'), default=sys.stdout, required=False)
@click.option('--exac-fields', '-e', type=str, default=None, help="csv of exac fields to include")
@click.option('--exac-field-default', '-d', type=(str, str), multiple=True, default=None, required=False, help='supply default value for exac fields')
def annotate(filename, output, exac_fields, exac_field_default):
	"""
Annotate entries in a vcf file\n
FILENAME required; path to an input vcf file\n
OUTPUT optional; path to output, will default to sys.stdout
"""
	if exac_fields is not None:
		exac_fields = exac_fields.split(',')
	else:
		exac_fields = []
	defs = dict()
	for item in exac_field_default:
		if item[0] in exac_defaults:
			raise click.BadParameter(
				"Attempted to define default for %s which is automatically supplied, existing value: %s new value: %s" % (item[0], exac_defaults[item[0]], item[1]))
		if item[0] in defs:
			raise click.BadParameter(
				"Defined default for %s twice, existing value: %s new value: %s" % (item[0], defs[item[0]], item[1]))
		if item[0] not in exac_fields:
			raise click.BadParameter(
				"Defined default for %s, but this field is not defined in --exac-fields" % (item[0]))
		defs[item[0]] = item[1]
	for field in exac_fields:
		if field in exac_vars:
			raise click.BadParameter(
				"Requested ExAC variable twice: %s" % field)
		if field in defs:
			add_exac_var(field, defs[field])
		else:
			add_exac_var(field)

	output_handle = output
	rdr = vcf.Reader(fsock=filename)
	# metakeys = list(rdr.metadata.keys())
	# print("META ENTRIES")
	# for key in metakeys:
	# 	print("%s\t%s" % (key, rdr.metadata[key]))
	# formatkeys = list(rdr.formats.keys())
	# print("FORMAT ENTRIES")
	# for key in formatkeys:
	# 	fmt = rdr.formats[key]
	# 	print("%s\t%s\t%s" % (key, fmt.type, fmt.desc))
	# infokeys = list(rdr.infos.keys())
	# print("INFO ENTRIES")
	# for key in infokeys:
	# 	inf = rdr.infos[key]
	# 	print("%s\t%s\t%s" % (key, inf.type, inf.desc))
	# entry = next(rdr)
	# read the vcf entries from the file
	keystrings = []  # type: List[str]
	records = []  # type: List[AnnotationRecord]
	#get records from vcf file
	for entry in rdr:
		annrec = AnnotationRecord.from_vcf_entry(entry)
		for ann in annrec:
			records.append(ann)
			keystrings.append(ann.varkey)
	#pull json data from ExAC
	json_data = ExAC.get_bulk_variant_data(keystrings)
	#use json data to fill values into records
	for i in range(len(keystrings)):
		keystr = keystrings[i]
		record = records[i]
		json_record = json_data[keystr]
		record.exac_data = ExAC_VariantData(exac_vars, exac_defaults, json_record)
	#print filled records
	exac_var_headers = []
	for var in exac_vars:
		exac_var_headers.append("ExAC_%s" % (var.replace('.', '_')))
	cols = ["chrom", "pos", "ref", "alt", "depth", "ref_supporting_reads", "ref_support_percent", "var_supporting_reads", "var_support_percent", "\t".join(exac_var_headers)]
	print("\t".join(cols), file=output_handle)
	for record in records:
		print(str(record), file=output_handle)

if __name__ == "__main__":
	annotate()
