#!/usr/bin/env python3

"""
Author: Bram van Wersch
University: Wageningen university
Date: 13/06/2019
"""

import tsv_analyser as tsva
import argparse

def get_arguments():
	"""
	Function using argparse to parse command line arguments.
	"""
	parser = argparse.ArgumentParser(description='Add dispersed duplications as insertions to the vcf file so vep takes these into account.')
	parser.add_argument("-vcf", help = "Input file in VCF format.", required = True)
	parser.add_argument("-o", "--output", help = "Loaction for the output file.", required = True)
	return parser.parse_args()

class AddDispercedDuplications:
	def __init__(self, vcf_file, out_file):
		self.vcf_tsv = tsva.TabSeperatedValueAnalyser2(vcf_file, "ID")
		self.add_disperced_insertions()
		self.add_to_header()
		self.vcf_tsv.vcf_to_file(out_file, "INFO")

	def add_disperced_insertions(self):
		ins_list = []
		for key in self.vcf_tsv.keys():
			if "SVTYPE=DUP:DISPERSED" in self.vcf_tsv[key]["INFO"]:
				ins_row = self.create_insertion_row(self.vcf_tsv[key])
				ins_list.append(ins_row)
		#making sure dictionary size is not changed while itterating.
		for row in ins_list:
			self.vcf_tsv[row["ID"]] = row

	def create_insertion_row(self, row):
		new_row = row.copy()
		#create a dictionary from the information in the info column for
		#easier navigation of that column
		info_dict = self.get_info_dict(row["INFO"])
		new_row["CHROM"] = info_dict["INSCHROM"]
		new_row["POS"] = info_dict["INSPOS"]
		new_row["ID"] = row["ID"] + ".i"
		new_row["ALT"] =  "<INS>" # cannot realy paste the whole insertion sequence here
		new_row["INFO"] = self.create_info_value(info_dict)
		return new_row
	
	def create_info_value(self, info_dict):
		info_list = []
		for key in info_dict:
			if key not in ["END","SVTYPE","INSCHROM","INSPOS"]:
				info_list.append("{}={}".format(key, info_dict[key]))
		info_list.append("SVTYPE=INS:DISPERSED")
		info_str = ";".join(info_list)
		return info_str
	
	def get_info_dict(self, info_row):
		info_dict = {}
		for value in info_row.split(";"):
			name, val = value.split("=")
			info_dict[name] = val
		return info_dict
		
	def add_to_header(self):
		#add description of insertion
		at_alt_header = False
		for x in range(len(self.vcf_tsv.header)):
			row = self.vcf_tsv.header[x]
			if row.startswith("##ALT="):
				at_alt_header = True
			if at_alt_header and not row.startswith("##ALT="):
				self.vcf_tsv.header.insert(x, "##ALT=<ID=INS:DISPERSED,Description=\"Insertion site of dispersed duplication\">")
				return	
		
#add iterator class for going over tab seperated files.
if __name__ == "__main__":
	args = get_arguments()
	AddDispercedDuplications(args.vcf, args.output)
