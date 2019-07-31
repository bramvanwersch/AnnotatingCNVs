#!/usr/bin/env python3

"""
Author: Bram van Wersch
University: Wageningen university
Date: 13/06/2019
"""

import sys


class VcfAnalyser:

    def __init__(self, file_loc, key_name=None):
        """
        Class that takes a vcf or vep file and stores it as dictionary of dictionaries for manipulation purposes.
        :param file_loc: The location of the vcf or vep file.
        :param key_name: Potential name of a column in the vcf or vep file to function as key for the saving of
        the inner dictionaries
        :param tsv_string: string containing the whole file.
        :param tsv_dict: the dictionary of dictionaries that holds all the data neccesairy.
        :param header: list that holds the header of the repsective file
        :param column_header: list that holds the column header of the repsective file.
        """
        self.file_loc = file_loc
        self.key_name = key_name
        self.tsv_string = self.get_tsv_text(file_loc)
        self.tsv_dict = {}
        self.header = []
        self.column_header = []
        #create the dictionary of dictionaries while innialising the class
        self.tsv_to_dict()
        
    def get_tsv_text(self, file_loc):
        """
        Function that opens a file and returns its contents
        :return: String that contains the whole file
        """
        t = open(file_loc)
        tsv_text = t.read()
        t.close()
        return tsv_text

    def tsv_to_dict(self):
        """
        Function that loops trough all lines in a vcf or vep file and seperates header, column header and the tab
        seperated rows. The rows get put into a dictionary of dictionaries with the outer dictionary containing key names
        as provided by self.key_name and the inner dictionary containing column names and their respective values.
        """
        tsv_rows = []
        for line in self.tsv_string.split("\n"):
            #header lines
            if line.startswith("##"):
                self.header.append(line)
            #column header line.
            elif line.startswith("#"):
                self.column_header = line[1:].split("\t")
            elif line:
                tsv_rows.append(line)
        for x in range(len(tsv_rows)):
            #seperates each row element seperated by tabs into a dictionary
            inner_dict = dict(zip(self.column_header, tsv_rows[x].split("\t")))
            #if not unique key is provided the outer dictionary is simply saved by a number.
            if self.key_name:
                self.tsv_dict[inner_dict[self.key_name]] = inner_dict
            else:
                self.tsv_dict[str(x)] = inner_dict

    def vep_to_file(self, out_file, info_key):
        """
        Function that takes a vep dictiory fo dictionaries and writes a sorted string to a file.
        :param out_file: Location and name of the file to write the output into
        :param info_key: key that contains the info column for furter sorting after location. This is relevant for
        dispersed duplications that can have the same location but are inserted in different locations.
        """
        tsv_text = self.to_string(info_key)
        t = open(out_file, "w")
        t.write(tsv_text)
        t.close()

    def vcf_to_file(self, out_file, info_key):
        """
        Function that takes a vcf dictiory fo dictionaries and writes a sorted string to a file. An extra step of
        creating a location value is needed for proper sorting of the values.
        :param out_file: Location and name of the file to write the output into
        :param info_key: key that contains the info column for furter sorting after location. This is relevant for
        dispersed duplications that can have the same location but are inserted in different locations.
        """
        #add a location key to the inner dictionary to sort values on
        for key, row in self.tsv_dict.items():
            loc_value = self.create_loc_value(row)
            self.tsv_dict[key]["Location"] = loc_value
        tsv_text = self.to_string(info_key)
        t = open(out_file, "w")
        t.write(tsv_text)
        t.close()

    def create_loc_value(self, row):
        """
        Function that creates a location value in the genome browser syntax. This is only relevant for vcf files that
        do not contain a nice location value to easily sort on.
        :param row: Dicitonary representing a row of the vcf file.
        :return: String containing a location value in the genome browser syntax.
        """
        end_val = [x for x in row["INFO"].split(";") if "END=" in x]
        try:
            end_loc = end_val[1]
        except IndexError:
            # in case of an insertion there is no end pos so just take the start
            end_loc = row["POS"]
        loc_value = "{}:{}-{}".format(row["CHROM"], row["POS"], *end_loc)
        return loc_value

    def to_string(self, info_key):
        """
        Function that creates a string from the self.tsv_dict value.
        :param info_key: key that contains the info column for furter sorting after location. This is relevant for
        dispersed duplications that can have the same location but are inserted in different locations.
        :return: A string containing the whole self.tsv_dict value together with the header and column header back into
        the csv or vep file it originally was. Sorted on location.
        """
        dict_list = [value for key, value in self.tsv_dict.items()]
        dict_list.sort(key=lambda x: self.sort_by(x["Location"], x[info_key]))
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
        Function for sorting creating tuples that the sort function of python can sort the rows on.
        :param loc: The value saved under the 'Location' key in the inner dict.
        :param info: The value saved undet he 'INFO' key in the inner dict.
        :return: A tuple containing True or false depending on if the chromosome was an integer or not. The chromonsome
        value, the start of the variant, the end of the variant, the info value associated with the variant.
        """
        chrom = loc.split(":")[0]
        if "-" in loc:
            start, end = loc.split(":")[1].split("-")
        # in case of insertions
        else:
            start = loc.split(":")[1]
            end = start
        # see if the chromosome can be made into an integer If this is the
        # case False is put at the start to make the numbered chromosomes come first
        try:
            return (False, int(chrom), int(start), int(end), info)
        except ValueError:
            return (True, chrom, int(start), int(end), info)

######## MAGIC METHODS #############
# these methods are here to make the class instance behave as a dictionary.

    def keys(self):
        return self.tsv_dict.keys()

    def __setitem__(self, key, value):
        self.tsv_dict[key] = value

    def __getitem__(self, key):
        return self.tsv_dict[key]

    def __delitem__(self, key):
        del self.tsv_dict[key]
