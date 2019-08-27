#!/usr/bin/env python3

"""
Author: Bram van Wersch
University: Wageningen university
Date: 13/06/2019
"""

import python_scripts.vcf_analyser as vcfa
import argparse


def get_arguments():
    """
    Function using argparse to parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Add dispersed duplications as insertions to the vcf file so vep takes these into account.')
    parser.add_argument("-vcf", help="Input file in VCF format.", required=True)
    parser.add_argument("-o", "--output", help="Loaction for the output file.", required=True)
    return parser.parse_args()


class AddDispercedDuplications:
    def __init__(self, vcf_file, out_file):
        """
        Class that takes a vcf file and adds dispersed duplications as insertions to the vcf file. This is done to make
        sure vep measures the impact of the dispersed insertions
        :param vcf_file: The input file in vcf format.
        :param out_file: The location and name of the output file.
        :param vcf_tsv: instance of VcfAnalyser class that holds the file as a dictionary of dictionaries and can be
        manipulated as such.
        """
        self.vcf_tsv = vcfa.VcfAnalyser(vcf_file, "ID")
        self.add_disperced_insertions()
        self.add_to_header()
        self.vcf_tsv.vcf_to_file(out_file, "INFO")

    def add_disperced_insertions(self):
        """
        Function that guides the process of adding the dispersed duplications
        """
        ins_list = []
        for key in self.vcf_tsv.keys():
            if "SVTYPE=DUP:DISPERSED" in self.vcf_tsv[key]["INFO"]:
                ins_row = self.create_insertion_row(self.vcf_tsv[key])
                ins_list.append(ins_row)
        # Creating a list of insertions and adding them after making sure the siuze of the dictionary does not change
        #while itterating over it.
        for row in ins_list:
            self.vcf_tsv[row["ID"]] = row

    def create_insertion_row(self, row):
        """
        Function that creates the row that has to be inserted based on the row it is creating the insertion from.
        :param row: Dictionary representing a dispersed insertion row that has has to create a new insertion.
        :return: Dictionary that has all values changed that should be changed for it to represnet the insertion the
        dispersed duplication would have caused.
        """
        # copy the row to prevent manipulating the existing row.
        new_row = row.copy()
        info_dict = self.get_info_dict(row["INFO"])
        new_row["CHROM"] = info_dict["INSCHROM"]
        new_row["POS"] = info_dict["INSPOS"]
        new_row["ID"] = row["ID"] + ".i"
        new_row["ALT"] = "<INS>"
        new_row["INFO"] = self.create_info_value(info_dict)
        return new_row


    def get_info_dict(self, info_string):
        """
        Fucntion that creates a dictionary out of the value saved under the 'INFO' key. This is done to make extracting
        the information from this value allot easier and more readable
        :param info_string: String that is saved under the 'INFO' key in the
        :return: A dictionary containing the values present in the info_string provided.
        """
        info_dict = {}
        for value in info_string.split(";"):
            name, val = value.split("=")
            info_dict[name] = val
        return info_dict

    def create_info_value(self, info_dict):
        """
        Function that creates the new INFO value for the dispersed insertion.
        :param info_dict: A dictionary form the INFO string of the original row.
        :return: A string that contains all relevant information that is part of the info string.
        """
        info_list = []
        for key in info_dict:
            if key not in ["END", "SVTYPE", "INSCHROM", "INSPOS"]:
                info_list.append("{}={}".format(key, info_dict[key]))
        info_list.append("SVTYPE=INS:DISPERSED")
        info_str = ";".join(info_list)
        return info_str

    def add_to_header(self):
        """
        Function that adds a line to the header giving a short explanation about the dispersed insertions.
        :return:
        """
        at_alt_header = False
        for x in range(len(self.vcf_tsv.header)):
            row = self.vcf_tsv.header[x]
            if row.startswith("##ALT="):
                at_alt_header = True
            if at_alt_header and not row.startswith("##ALT="):
                self.vcf_tsv.header.insert(x,"##ALT=<ID=INS:DISPERSED,Description=\"Insertion site of dispersed duplication\">")
                return

if __name__ == "__main__":
    args = get_arguments()
    AddDispercedDuplications(args.vcf, args.output)
