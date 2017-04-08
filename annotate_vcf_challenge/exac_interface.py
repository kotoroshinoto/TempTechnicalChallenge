from vcf.model import _Record as vcfrecord
import requests
import json
from typing import List
import annotate_vcf_challenge.json_func as json_func

consequence_severity = {
	'intergenic_variant': 0,
	'feature_truncation': 1,
	'regulatory_region_variant': 2,
	'feature_elongation': 3,
	'regulatory_region_amplification': 4,
	'regulatory_region_ablation': 5,
	'TF_binding_site_variant': 6,
	'TFBS_amplification': 7,
	'TFBS_ablation': 8,
	'downstream_gene_variant': 9,
	'upstream_gene_variant': 10,
	'non_coding_transcript_variant': 11,
	'NMD_transcript_variant': 12,
	'intron_variant': 13,
	'non_coding_transcript_exon_variant': 14,
	'3_prime_UTR_variant': 15,
	'5_prime_UTR_variant': 16,
	'mature_miRNA_variant': 17,
	'coding_sequence_variant': 18,
	'synonymous_variant': 19,
	'stop_retained_variant': 20,
	'incomplete_terminal_codon_variant': 21,
	'splice_region_variant': 22,
	'protein_altering_variant': 23,
	'missense_variant': 24,
	'inframe_deletion': 25,
	'inframe_insertion': 26,
	'transcript_amplification': 27,
	'start_lost': 28,
	'stop_lost': 29,
	'frameshift_variant': 30,
	'stop_gained': 31,
	'splice_donor_variant': 32,
	'splice_acceptor_variant': 33,
	'transcript_ablation': 34
}


class ExAC:
	@staticmethod
	def split_keystring(keystring:str):
		splitstr = keystring.split('-')
		chrom = splitstr[0]
		pos = int(splitstr[1])
		ref = splitstr[2]
		alt = splitstr[3]
		return chrom, pos, ref, alt

	@staticmethod
	def get_most_severe_conseqeuence(variant, json_data):
		json_entry = json_data[variant]
		vep_annotations = json_func.expand_json_var(json_entry, "variant.vep_annotations")
		if vep_annotations is None:
			return None
		worst = None
		worst_val = -1
		for annotation in vep_annotations:
			cons = annotation['major_consequence']
			if cons not in consequence_severity:
				raise ValueError("consequence not recognized: '%s'" % cons)
			cons_val = consequence_severity[cons]
			if cons_val > worst_val:
				worst = cons
				worst_val = cons_val
		return worst

	@staticmethod
	def normalize_keystring(keystr: str):
		chrom, pos, ref, alt = ExAC.split_keystring(keystr)
		num_trunc = 0
		for i in range(1, min(len(ref), len(alt))):
			if ref[len(ref) - i] == alt[len(alt)-i]:
				num_trunc += 1
				continue
			else:
				break
		if num_trunc > 0:
			ref = ref[:-num_trunc]
			alt = alt[:-num_trunc]
		num_trunc = 0
		for i in range(0, min(len(ref)-1, len(alt)-1)):
			if ref[i] == alt[i]:
				num_trunc += 1
			else:
				break
		if num_trunc > 0:
			ref = ref[num_trunc:]
			alt = alt[num_trunc:]
			pos += num_trunc
		return "%s-%d-%s-%s" % (chrom, pos, ref, alt)

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

	class VariantData:
		def __init__(self, varlist, defaults, json_data):
			self.vals = []
			#for each requested ExAC variable
			for var in varlist:
				#get required value from json data
				var_val = json_func.expand_json_var(json_data, var)
				#if the value is missing from json data and there is a default defined, use that
				if (var_val is None) and (var in defaults):
					self.vals.append(defaults[var])
				else:
					if isinstance(var_val, list):
						self.vals.append(",".join(str(x) for x in var_val))
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
