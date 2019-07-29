#!/usr/bin/env python3

"""
Author: Bram van Wersch
University: Wageningen university
Date: 13/06/2019
"""


import sys
import tsv_analyser as tsva
import argparse

def get_arguments():
	"""
	Function using argparse to parse command line arguments.
	"""
	parser = argparse.ArgumentParser(description='Add an extra information column to a generated VEP file.')
	parser.add_argument("--vcf", help = "Input file in VCF format.", required = True)
	parser.add_argument("--vep", help = "VEP output file.", required = True)
	parser.add_argument("-o", "--output", help = "Loaction for the output file.", required = True)
	parser.add_argument("--vepIDName", help = "Name of VEP identifier column. Default = 'Uploaded_variation'.", required = False, default = "Uploaded_variation")
	parser.add_argument("--vcfIDName", help = "Name of VCF identifier column. Default = 'ID'.", required = False, default = "ID")
	parser.add_argument("--vcfColumnNames", nargs = "+", help = "Names of VCF column(s) to be added to VEP file. Default = 'INFO'.", required = False, default = ["INFO"])
	return parser.parse_args()

class AddInformationVep:
	"""
	Class that will when called take the provided vep and vcf file and add
	columns from the vcf file and add them to the vep file.
	vcf_in -- string location of the vcf file to be used
	vep_in -- string location of the vep file to be used
	out_file -- string location of the output file to be placed.
	vep_ID_name -- string name of the column that contains the copy number variant
	identifiers for the vep file
	vcf_ID_name -- string name of the column that contains the copy number variant
	identifiers for the vcf file
	vcf_column_names -- list containing names of all columns that have to be 
	copied from the vcf file to the vep file.
	"""
	def __init__(self, vcf_in, vep_in, out_file, vep_ID_name, vcf_ID_name, vcf_column_names):
		self.vcf_in = vcf_in
		self.vep_in = vep_in
		self.out_file = out_file
		self.vep_ID_name = vep_ID_name
		self.vcf_ID_name = vcf_ID_name
		self.vcf_column_names = vcf_column_names
		
	
	def protocol(self):
		"""
		Function that runs the protocol for this class.
		1. Create 2 objects of TabSeperatedValueAnalyser one for each input
		file.
		2. Add the columns to the TabSeperatedValueAnalyser object of the
		vep file.
		3. Add to the vep header to make it contain all needed information.
		4. Write the output to the requested output file.
		"""
		#1.
		self.vcf_tsv = tsva.TabSeperatedValueAnalyser2(self.vcf_in, self.vcf_ID_name)
		self.vep_tsv = tsva.TabSeperatedValueAnalyser2(self.vep_in)
		#2.
		self.add_information()
		#3.
		self.add_to_header()
		#4.
		self.vep_tsv.vep_to_file(self.out_file, "Extra")
	
	def add_information(self):
		"""
		Function that extends dictionary entries with an INFO column and 
		adds some extra information to certain variants that lack that
		information to the self.vep_tsv dictionarie. 
		"""
		#add the information that is going to be added to the column header
		col_names = [name for name in self.vcf_column_names]
		self.vep_tsv.column_header += col_names
		for row_key in self.vep_tsv.keys():
			vep_ID = self.vep_tsv[row_key][self.vep_ID_name]
			self.check_insertion_info(self.vep_tsv[row_key])
			self.check_start_codon(self.vep_tsv[row_key])
			for name in self.vcf_column_names:
				try:
					vcf_info = self.vcf_tsv[vep_ID][name]	
					self.vep_tsv[row_key][name] = vcf_info
				except KeyError:
					#better error return
					#vep ID that does not point to a correct line but empty line or incomplete line.
					#also remove the line.
					print(vep_ID)
					
	def check_insertion_info(self, row_dict):
		"""
		Function that checks if a certain variant is an insertion and if it causes a 
		frame shift. If this is the case additional information is added to make sure
		it is annotated that the insertion is a frame shift.
		row_dict -- dictionary entrie representing one row in the vep file.
		TODO: somehow check for introduction of start codons
		"""
		if row_dict["Allele"] == "insertion" and "coding_sequence_variant" in row_dict["Consequence"]:
			insertion = self.vcf_tsv[row_dict[self.vep_ID_name]]["ALT"]
			if insertion != "<INS>" and (len(insertion) - 1) % 3 != 0:
				row_dict["Consequence"] += ",frameshift_variant"
				## replace the consequence with high. 
				row_dict["Extra"] = row_dict["Extra"].replace("MODIFIER", "HIGH").replace("LOW", "HIGH").replace("MODERATE", "HIGH")
			elif insertion != "<INS>" and (len(insertion) - 1) % 3 == 0:
				row_dict["Consequence"] += ",inframe_insertion"
				## replace the consequence with moderate --> this is what the vep documentation states.. 
				row_dict["Extra"] = row_dict["Extra"].replace("MODIFIER", "MODERATE").replace("LOW", "MODERATE")

				
				
	def check_start_codon(self, row_dict):
		"""
		Function that checks if there should have been a start_lost annotated. Or removed if the 
		allele is a duplication. This is to fix the broken start_lost consequence annotated by VEP.
		"""
		if row_dict["Allele"] == "deletion" and "5_prime_UTR_variant" in row_dict["Consequence"] and\
			("coding_sequence_variant" or "frameshift_variant") in row_dict["Consequence"] and "start_lost" not in row_dict["Consequence"]:
			row_dict["Consequence"] += ",start_lost"
		elif row_dict["Allele"] == "duplication" and "start_lost" in row_dict["Consequence"]:
			row_dict["Consequence"].replace(",start_lost", "")
			
		
	def add_to_header(self):
		"""
		Function that adds info from the vcf_header to the vep_header depending
		on the columns added and adds the vep_column info together with values 
		for the additional columns to the vep_header.
		"""
		for line in self.vcf_tsv.header:
			for name in self.vcf_column_names:
				if line.startswith("##" + name):
					self.vep_tsv.header.append(line)


if __name__ == "__main__":
	#add check to see if correct file types
	#add check to see if file location already is in use
	args = get_arguments()
	AddInformationVep(args.vcf, args.vep, args.output, args.vepIDName, args.vcfIDName, args.vcfColumnNames).protocol()
	
	
