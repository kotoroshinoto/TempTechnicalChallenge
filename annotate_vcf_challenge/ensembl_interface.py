import requests
import json
from typing import List
import annotate_vcf_challenge.json_func as json_func
from annotate_vcf_challenge.exac_interface import ExAC


class Ensembl:
	base_url = "http://grch37.rest.ensembl.org"

	@staticmethod
	def convert_exac_keystring_to_url(keystr:str):
		chrom, pos, ref, alt = ExAC.split_keystring(keystr)
		endpos = pos + len(ref) - 1
		return "%s:%d-%d/%s" % (chrom, pos, endpos, alt)

	@staticmethod
	def convert_exac_keystring_for_post(keystr: str):
		chrom, pos, ref, alt = ExAC.split_keystring(keystr)
		return "%s %d . %s %s" % (chrom, pos, ref, alt)

	@staticmethod
	def get_vep_variant_url(urlkeystring):
		url_suffix="vep/homo_sapiens/region"
		return "%s/%s/%s" % (Ensembl.base_url, url_suffix, urlkeystring)

	@staticmethod
	def get_vep_variant_bulk_url():
		url_suffix = "vep/homo_sapiens/region"
		return "%s/%s" % (Ensembl.base_url, url_suffix)

	@staticmethod
	def get_variant_data(urlkeystring):
		headers = {'Accept': 'application/json'}
		url = Ensembl.get_vep_variant_url(urlkeystring)
		return requests.get(url, headers=headers)

	@staticmethod
	def get_bulk_variant_data(keystr_list: 'List[str]', verbose=False):
		headers ={"Content-Type": "application/json", 'Accept': 'application/json'}
		url = Ensembl.get_vep_variant_bulk_url()
		#can only handle 300 at once
		result = []
		buff = []
		numrec = len(keystr_list)
		i = 0
		for i in range(numrec):
			#we have 300 in buff, perform query and reset buff
			if i % 300 == 0:
				if len(buff) > 0:
					print("record %d of %d" % (i, numrec)) if verbose else None
					data = json.dumps({'variants': buff})
					r = requests.post(url, headers=headers, data=data)
					j = r.json()
					if not isinstance(j, list):
						raise RuntimeError("Did not get a list from json request! There may have been an error")
					for element in j:
						result.append(element)
				buff = []
			buff.append(keystr_list[i])
		#perform last query
		if len(buff) > 0:
			print("record %d of %d" % (i + 1, numrec)) if verbose else None
			data = json.dumps({'variants': buff})
			r = requests.post(url, headers=headers, data=data)
			j = r.json()
			if not isinstance(j, list):
				raise RuntimeError("Did not get a list from json request! There may have been an error")
			for element in j:
				result.append(element)
		return result

