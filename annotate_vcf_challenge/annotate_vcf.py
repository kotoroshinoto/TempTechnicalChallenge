#!/usr/bin/env python3
import sys
from typing import List
import shutil
import click
import vcf
from vcf.model import _Record as vcfrecord

from annotate_vcf_challenge.exac_interface import ExAC

output_handle = None
exac_vars = list()
exac_defaults = dict()


def add_exac_var(varstr, defval=None):
	exac_vars.append(varstr)
	if defval is not None:
		exac_defaults[varstr] = defval

# TODO decide if any additional information available from EXAC is relevant enough to add by default
#by default we want the allele frequency
add_exac_var('variant.allele_freq', 0.0)


class AnnotationRecord:
	def __init__(self, varstring: str):
		self.varkey = varstring
		self.depth = 0
		self.ref_support = 0
		self.ref_support_percent = 0.0
		self.var_support = 0
		self.var_support_percent = 0.0
		self.exac_data = None

	def __str__(self):
		cols = [self.varkey.replace('-',"\t"), self.depth,
		        self.ref_support, self.ref_support_percent,
		        self.var_support, self.var_support_percent]
		if self.exac_data is not None:
			cols.append(self.exac_data)
		return "\t".join(str(x) for x in cols)

	@classmethod
	def from_vcf_entry(cls, entry: vcfrecord):
		#keystrings are of the format ExAC expects for variants: chrom-pos-ref-alt
		keystr_list = ExAC.make_keystring_for_vcfentry(entry)
		records = []
		#create one entry per alt, otherwise we're treating multiple variant definitions as if they were a single variant
		for i in range(len(entry.ALT)):
			keystr = ExAC.normalize_keystring(keystr_list[i])
			record = cls(keystr)
			# print(keystr, file=output_handle)
			# print(entry.INFO, file=output_handle)
			# print(entry.samples, file=output_handle)
			#get read depth for variant location:
			record.depth = entry.INFO['DP']
			#get read observation count
			record.ref_support = entry.INFO['RO']
			#calculate read percentage contributing to ref
			record.ref_support_percent = float(record.ref_support) / float(record.depth)
			#get alt observation count
			record.var_support = entry.INFO['AO'][i]
			# calculate read percentage contributing to alt
			record.var_support_percent = float(record.var_support) / float(record.depth)
			# TODO determine variant type
			#add record to list
			records.append(record)
		#return record list
		return records


def handle_field_args(exac_fields, exac_field_default):
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


def read_vcf(file_handle):
	#use vcf reader from pyvcf
	rdr = vcf.Reader(fsock=file_handle)
	# read the vcf entries from the file
	keystrings = []  # type: List[str]
	records = []  # type: List[AnnotationRecord]
	# get records from vcf file
	for entry in rdr:
		records_from_vcf_entry = AnnotationRecord.from_vcf_entry(entry)
		for record in records_from_vcf_entry:
			records.append(record)
			#keystrings are of the format ExAC expects for variants: chrom-pos-ref-alt
			keystrings.append(record.varkey)
	return keystrings, records


def write_result(records):
	#print filled records
	exac_var_headers = []
	for var in exac_vars:
		exac_var_headers.append("ExAC_%s" % (var.replace('.', '_')))
	cols = ["chrom", "pos", "ref", "alt", "depth", "ref_supporting_reads", "ref_support_percent", "var_supporting_reads", "var_support_percent", "\t".join(exac_var_headers)]
	print("\t".join(cols), file=output_handle)
	for record in records:
		print(str(record), file=output_handle)


def update_records_with_exac_data(keystrings, records, json_data):
	for i in range(len(keystrings)):
		keystr = keystrings[i]
		record = records[i]
		json_record = json_data[keystr]
		record.exac_data = ExAC.VariantData(exac_vars, exac_defaults, json_record)


@click.command(context_settings=dict(max_content_width=shutil.get_terminal_size().columns))
@click.argument('filename', type=click.File('r'))
@click.argument('output', type=click.File('w'), default=sys.stdout, required=False)
@click.option('--exac-fields', '-e', type=str, default=None, help="csv of additional exac fields to include, formatted: parent.child.grandchild, as if from: http://exac.hms.harvard.edu/rest/variant")
@click.option('--exac-field-default', '-d', type=(str, str), multiple=True, default=None, required=False, help='supply default value for exac field:  <FIELD VALUE>')
def annotate(filename, output, exac_fields, exac_field_default):
	"""
Annotate entries in a vcf file\n
FILENAME required; path to an input vcf file\n
OUTPUT optional; path to output, will default to sys.stdout
"""
	global output_handle
	output_handle = output

	#add any additional requested ExAC fields
	handle_field_args(exac_fields, exac_field_default)

	#read in vcf contents
	keystrings, records = read_vcf(filename)

	#pull json data from ExAC
	json_data = ExAC.get_bulk_variant_data(keystrings)

	#use json data to fill values into records
	update_records_with_exac_data(keystrings, records, json_data)

	#output records
	write_result(records)

if __name__ == "__main__":
	annotate()
