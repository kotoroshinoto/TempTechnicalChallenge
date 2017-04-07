from vcf.model import _Record as vcfrecord
import requests
import json
from typing import List


class ExAC:
	# turn vcf fields into an identifier string for the variant
	@staticmethod
	def make_keystring_for_vcfentry(entry: vcfrecord):
		return ["%s-%s-%s-%s" % (entry.CHROM, entry.POS, entry.REF, x) for x in entry.ALT]

	exac_base_url = "http://exac.hms.harvard.edu"

	@staticmethod
	def variant_url(keystring, specifier=None):
		if specifier is None:
			return "%s%s%s" % (ExAC.exac_base_url, "/rest/variant/", keystring)
		else:
			return "%s%s%s" % (ExAC.exac_base_url, ("/rest/variant/%s/" % specifier), keystring)

	@staticmethod
	def bulk_variant_url(specifier=None):
		if specifier is None:
			return "%s%s" % (ExAC.exac_base_url, "/rest/bulk/variant")
		else:
			return "%s%s" % (ExAC.exac_base_url, ("/rest/bulk/variant/%s" % specifier))

	@staticmethod
	def expand_json_var(json_data, var: str):
		#string of format parent.child.grandchild split into keys for json dict
		levels = var.split(".")
		curr = json_data
		#index into successive dicts
		for lvl in levels:
			if lvl in curr:
				curr = curr[lvl]
			else:
				#if target key doesn't exist, return None, caller will handle condition
				return None
		#after all levels are traversed, we have the value we need
		return curr

	class VariantData:
		def __init__(self, varlist, defaults, json_data):
			self.vals = []
			#for each requested ExAC variable
			for var in varlist:
				#get required value from json data
				var_val = ExAC.expand_json_var(json_data, var)
				#if the value is missing from json data and there is a default defined, use that
				if (var_val is None) and (var in defaults):
					self.vals.append(defaults[var])
				else:
					self.vals.append(var_val)

		def __str__(self):
			return "\t".join(str(x) for x in self.vals)

	@staticmethod
	def get_variant_data(keystring):
		return requests.get(ExAC.variant_url(keystring)).json()

	@staticmethod
	def get_bulk_variant_data(keystr_list: 'List[str]'):
		return requests.post(ExAC.bulk_variant_url(), data=json.dumps(keystr_list)).json()
