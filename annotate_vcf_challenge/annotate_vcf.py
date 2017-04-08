#!/usr/bin/env python3
import sys
from typing import List
import shutil
import click
import vcf
from vcf.model import _Record as vcfrecord
from annotate_vcf_challenge.exac_interface import ExAC
from annotate_vcf_challenge.ensembl_interface import Ensembl
import json
import os.path
from pprint import pprint
OUTPUT = None
VERBOSE = None
CACHE = None
REFRESH = None
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
		self.type = None
		self.most_severe_consequence = None
		self.exac_data = None

	@staticmethod
	def get_column_headings():
		exac_var_headers = []
		for var in exac_vars:
			exac_var_headers.append("ExAC_%s" % (var.replace('.', '_')))
		return "\t".join(["chrom", "pos", "ref", "alt", "depth", "ref_supporting_reads", "ref_support_percent",
		                  "var_supporting_reads", "var_support_percent", "variant_type", "most_severe_consequence",
		                  "\t".join(exac_var_headers)])

	def __str__(self):
		cols = [self.varkey.replace('-',"\t"), self.depth,
		        self.ref_support, self.ref_support_percent,
		        self.var_support, self.var_support_percent, self.type]
		if self.most_severe_consequence is None:
			cols.append('-')
		else:
			cols.append(self.most_severe_consequence)
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
			#need to adjust keystring to remove excess matching text, allowing exac to match properly
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
			#get INFO type field
			record.type = entry.INFO['TYPE'][i]
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
	sys.stderr.write("Reading VCF ...") if VERBOSE else None
	#use vcf reader from pyvcf
	rdr = vcf.Reader(fsock=file_handle)
	# read the vcf entries from the file
	keystrings = []  # type: List[str]
	ensembl_keystrings = []  # type: List[str]
	records = []  # type: List[AnnotationRecord]
	# get records from vcf file
	for entry in rdr:
		records_from_vcf_entry = AnnotationRecord.from_vcf_entry(entry)
		for record in records_from_vcf_entry:
			records.append(record)
			#keystrings are of the format ExAC expects for variants: chrom-pos-ref-alt
			keystrings.append(record.varkey)
			ensembl_keystrings.append(Ensembl.convert_exac_keystring_for_post(record.varkey))
	print(" done", file=sys.stderr) if VERBOSE else None
	return keystrings, ensembl_keystrings, records


def write_result(records):
	#print filled records
	sys.stderr.write("Writing table to file: '%s' ..." % OUTPUT.name) if VERBOSE else None
	cols = AnnotationRecord.get_column_headings()
	print(cols, file=OUTPUT)
	for record in records:
		print(str(record), file=OUTPUT)
	print(" done", file=sys.stderr) if VERBOSE else None


def update_records_with_exac_data(keystrings, records, json_data):
	sys.stderr.write("integrating ExAC data into records ...") if VERBOSE else None
	for i in range(len(keystrings)):
		keystr = keystrings[i]
		record = records[i]  # type: AnnotationRecord
		json_record = json_data[keystr]
		record.exac_data = ExAC.VariantData(exac_vars, exac_defaults, json_record)
		# record.most_severe_consequence = ExAC.get_most_severe_conseqeuence(keystr, json_data)
	print(" done", file=sys.stderr) if VERBOSE else None


def update_records_with_ensembl_data(keystrings, records, json_data):
	sys.stderr.write("integrating Ensembl VEP data into records ...") if VERBOSE else None
	for i in range(len(keystrings)):
		keystr = keystrings[i]  # type: str
		record = records[i]  # type: AnnotationRecord
		json_record = json_data[i]
		if json_record['input'] != keystr:
			raise RuntimeError("record input doesn't match keystring")
		record.most_severe_consequence = json_record['most_severe_consequence']
	print(" done", file=sys.stderr) if VERBOSE else None


def obtain_exac_data(keystrings: 'List[str]'):
	exac_cache_fname = 'last_exac.json'
	if CACHE and os.path.isfile(exac_cache_fname) and (not REFRESH):
		sys.stderr.write("Reading ExAC data from cache ...") if VERBOSE else None
		infile = open(exac_cache_fname, 'r')
		json_data = json.load(infile)
		infile.close()
	else:
		sys.stderr.write("Getting ExAC data from server ...") if VERBOSE else None
		json_data = ExAC.get_bulk_variant_data(keystrings)
		print(" done", file=sys.stderr) if VERBOSE else None
		if CACHE:
			sys.stderr.write("Writng ExAC data to cache ...") if VERBOSE else None
			outfile = open(exac_cache_fname, 'w')
			json.dump(json_data, outfile)
			outfile.close()
	print(" done", file=sys.stderr) if VERBOSE else None
	return json_data


def obtain_ensembl_data(ensembl_keystrings: 'List[str]'):
	ensembl_cache_fname = 'last_ensembl_vep.json'
	if CACHE and os.path.isfile(ensembl_cache_fname) and (not REFRESH):
		sys.stderr.write("Reading Ensembl VEP data from cache ...") if VERBOSE else None
		infile = open(ensembl_cache_fname, 'r')
		ensembl_json_data = json.load(infile)
		infile.close()
	else:
		print("Getting Ensembl VEP data from server ...", file=sys.stderr) if VERBOSE else None
		ensembl_json_data = (Ensembl.get_bulk_variant_data(ensembl_keystrings, verbose=VERBOSE))
		print("... done", file=sys.stderr) if VERBOSE else None
		if CACHE:
			sys.stderr.write("Writing Ensembl VEP data to cache ...") if VERBOSE else None
			outfile = open(ensembl_cache_fname, 'w')
			json.dump(ensembl_json_data, outfile)
			outfile.close()
	print(" done", file=sys.stderr) if VERBOSE else None
	return ensembl_json_data


def handle_globvars(output, verbose, cache, refresh):
	global OUTPUT
	global VERBOSE
	global CACHE
	global REFRESH
	OUTPUT = output
	VERBOSE = verbose
	CACHE = cache
	REFRESH = refresh


@click.command(context_settings=dict(max_content_width=shutil.get_terminal_size().columns))
@click.argument('filename', type=click.File('r'))
@click.argument('output', type=click.File('w'), default=sys.stdout, required=False)
@click.option('--exac-fields', '-e', type=str, default=None, help="csv of additional exac fields to include, formatted: parent.child.grandchild, as if from: http://exac.hms.harvard.edu/rest/variant")
@click.option('--exac-field-default', '-d', type=(str, str), multiple=True, default=None, required=False, help='supply default value for exac field:  <FIELD VALUE>')
@click.option('--cache', '-c', default=False, is_flag=True, help="if cache doesn't exist yet, save server response to file, if cache exists, read it instead")
@click.option('--verbose', '-v', default=False,  is_flag=True, help="verbose output to stderr")
@click.option('--refresh', '-r', default=False, is_flag=True, help="force a refresh of cache")
def annotate(filename, output, exac_fields, exac_field_default, cache, verbose, refresh):
	"""
Annotate entries in a vcf file\n
FILENAME required; path to an input vcf file\n
OUTPUT optional; path to output, will default to sys.stdout
"""
	#handle output filehandle and global flags
	handle_globvars(output, verbose, cache, refresh)

	#add any additional requested ExAC fields
	handle_field_args(exac_fields, exac_field_default)

	#read in vcf contents
	keystrings, ensembl_keystrings, records = read_vcf(filename)

	#pull json data from ExAC
	json_data = obtain_exac_data(keystrings)

	#use json data to fill values into records
	update_records_with_exac_data(keystrings, records, json_data)

	#get variant prediction data from ensembl
	ensembl_json_data = obtain_ensembl_data(ensembl_keystrings)

	#use json data to supply variant effect predictions
	update_records_with_ensembl_data(ensembl_keystrings, records, ensembl_json_data)

	#output records
	write_result(records)

if __name__ == "__main__":
	annotate()
