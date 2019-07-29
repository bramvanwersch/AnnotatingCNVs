#!/usr/bin/env python3

"""
Author: Bram van Wersch
University: Wageningen university
Date: 13/06/2019
"""
import sys

class TabSeperatedValueAnalyser2:
	"""
	Class that takes a file location and potentialy a key name from an
	information column to save the rows in the file by.
	file_loc -- Location of the tab seperated file.
	key_name -- name of a key to save each row of the tab seperated values by
	this name refers to a value in that row and can only be used if it is unique.
	"""
	def __init__(self, file_loc, key_name = None):
		self.file_loc = file_loc
		self.tsv_string = self.get_tsv_text(file_loc)
		self.key_name = key_name
		self.tsv_dict = {}
		self.header = []
		self.column_header = []
		self.tsv_to_dict()
				
	def get_tsv_text(self, file_loc):
		"""
		Function that opens a file and returns its contents
		return -- string containing the full contents of the file.
		"""
		t = open(file_loc)
		tsv_text = t.read()
		t.close()
		return tsv_text
	
	def tsv_to_dict(self):
		"""
		Function that goes trough all lines in the file and seperated the 
		header the column header and the tab seperated rows. The rows are
		then saved in a dictionary. The rows itself are represented as dictionaries
		where each row value is a dictionary entry. The names for the values
		are derived from the column_header.
		"""
		tsv_rows = []
		for line in self.tsv_string.split("\n"):
			if line.startswith("##"):
				self.header.append(line)
			elif line.startswith("#"):
				self.column_header = line[1:].split("\t")
			elif line:
				tsv_rows.append(line)	
		for x in range(len(tsv_rows)):
			inner_dict = dict(zip(self.column_header, tsv_rows[x].split("\t")))
			## takes the value of a key from the inner dict to save the row by.
			##this significantly speeds up things because you can directly ask for
			##a certain ID instead of looping over keys.
			if self.key_name:
				self.tsv_dict[inner_dict[self.key_name]] = inner_dict
			else:
				self.tsv_dict[str(x)] = inner_dict
	
	def vep_to_file(self, out_file, info_key):
		"""
		Function that writes a string to the specified out_file location
		is specified. Otherwise puts it in the same directory as the vep 
		file with an extra_info extension to the file name.
		"""
		tsv_text = self.to_string(info_key)
		t = open(out_file, "w+")
		t.write(tsv_text)
		t.close()
		
	def vcf_to_file(self, out_file, info_key):
		"""
		"""
		# for properly sorting them on location a location value has to be 
		#constructed. At this moment the locations are in three seperate columns
		for key,row in self.tsv_dict.items():
			loc_value = self.create_loc_value(row)
			self.tsv_dict[key]["Location"] = loc_value
		tsv_text = self.to_string(info_key)
		t = open(out_file, "w+")
		t.write(tsv_text)
		t.close()
		
	def create_loc_value(self, row):
		"""
		"""
		end_val = [x for x in row["INFO"].split(";") if "END=" in x]
		try:
			end_loc = end_val[1]
		except IndexError:
			#in case of an insertion there is no end pos so just take the start
			end_loc = row["POS"]
		loc_value = "{}:{}-{}".format(row["CHROM"], row["POS"], *end_loc)
		return loc_value
	
	def to_string(self, info_key):
		"""
		"""
		#sawp ID key with loc key
		dict_list = [value for key, value in self.tsv_dict.items()]
		dict_list.sort(key = lambda x: self.sort_by(x["Location"], x[info_key]))
		dict_str = ""
		for dictionary in dict_list:
			for name in self.column_header:
				dict_str += dictionary[name] + "\t"
			dict_str += "\n"
		self.header = "\n".join(self.header) + "\n"
		self.column_header = "#" + "\t".join(self.column_header) + "\n"
		return self.header + self.column_header + dict_str
			
	def sort_by(self, loc, info):
		"""
		"""
		chrom = loc.split(":")[0]
		if "-" in loc:
			start, end = loc.split(":")[1].split("-")
		# in case of insertions
		else:
			start = loc.split(":")[1]
			end = start
		# see if the chromosome can be made into an integer If this is the
		#case False is put at the start to make the numbered chromosomes come first
		try:
			return (False, int(chrom), int(start), int(end), info)
		except ValueError:
			return (True, chrom, int(start), int(end), info)
			
	def keys(self):
		return self.tsv_dict.keys()
	
	def __setitem__(self, key, value):
		self.tsv_dict[key] = value
		
	def __getitem__(self, key):
		return self.tsv_dict[key]
	
	def __delitem__(self, key):
		del self.tsv_dict[key]

if __name__ == "__main__":
	arg = sys.argv[1]
	TabSeperatedValueAnalyser2(arg)
