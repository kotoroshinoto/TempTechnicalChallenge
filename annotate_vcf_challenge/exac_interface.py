from vcf.model import _Record as vcfrecord
import requests
import json
from typing import List


# turn vcf fields into an identifier string for the variant
def make_keystring_for_vcfentry(entry: vcfrecord):
	return ["%s-%s-%s-%s" % (entry.CHROM, entry.POS, entry.REF, x) for x in entry.ALT]


class ExacUrlConstructor:
	base_url = "http://exac.hms.harvard.edu"

	@staticmethod
	def variant(keystring, specifier=None):
		if specifier is None:
			return "%s%s%s" % (ExacUrlConstructor.base_url, "/rest/variant/", keystring)
		else:
			return "%s%s%s" % (ExacUrlConstructor.base_url, ("/rest/variant/%s/" % specifier), keystring)

	@staticmethod
	def bulk_variant(specifier=None):
		if specifier is None:
			return "%s%s" % (ExacUrlConstructor.base_url, "/rest/bulk/variant")
		else:
			return "%s%s" % (ExacUrlConstructor.base_url, ("/rest/bulk/variant/%s" % specifier))


def expand_json_var(json_data, var: str):
	levels = var.split(".")
	curr = json_data
	for lvl in levels:
		# print("attempting to expand using: %s" % lvl, file=sys.stderr)
		if lvl in curr:
			curr = curr[lvl]
		else:
			# print("%s not found in json dict" % lvl, file=sys.stderr)
			return None
	# print("ACTUAL VALUE RETURNED FROM EXPAND: %s" % curr, file=sys.stderr)
	return curr


class ExAC_VariantData:
	def __init__(self, varlist, defaults, json_data):
		self.vals = []
		for var in varlist:
			var_val = expand_json_var(json_data, var)
			if (var_val is None) and (var in defaults):
				self.vals.append(defaults[var])
			else:
				self.vals.append(var_val)

	def __str__(self):
		return "\t".join(str(x) for x in self.vals)


class ExAC:
	@staticmethod
	def get_variant_frequency(keystring):
		r = requests.get(ExacUrlConstructor.variant(keystring, specifier="variant"))
		j = r.json()
		if "allele_freq" in j:
			# pprint(j['allele_freq'])
			return j['allele_freq']
		else:
			return 0.0

	@staticmethod
	def get_bulk_variant_data(keystr_list: 'List[str]'):
		# print(keystr_list)
		r = requests.post(ExacUrlConstructor.bulk_variant(), data=json.dumps(keystr_list))
		# print(r.text)
		j = r.json()
		return j
		# for keystr in keystr_list:
		# 	# pprint(j[keystr])
		# 	if 'allele_freq' in j[keystr]:
		# 		freqs.append(j[keystr]['allele_freq'])
		# 	else:
		# 		freqs.append(0.0)
		# return freqs

