#!/usr/bin/env python3

"""
Author: Bram van Wersch
University: Wageningen university
Date: 13/06/2019
"""

import sys
import python_scripts.vcf_analyser as vcfa
import argparse

def get_arguments():
    """
    Function using argparse to parse command line arguments.
    """
    parser = argparse.ArgumentParser(description='Add an extra information column to a generated VEP file.')
    parser.add_argument("--vcf", help="Input file in VCF format.", required=True)
    parser.add_argument("--vep", help="VEP output file.", required=True)
    parser.add_argument("-o", "--output", help="Loaction for the output file.", required=True)
    return parser.parse_args()

class AddInformationVep:
    def __init__(self, vcf_in, vep_in, out_file):
        """
        Class that checks the vep file for problems that where found scanning the output. These problems are: start
        codons not being annotated correctly, duplications causing start_lost variants and insertions not causing a
        frame shift or inframe insertion. Besides that the INFO column from the vcf file is added to the vep file
        because it contains unique information that would get lost otherwise.
        :param vcf_in: location of the vcf file that has to have its 'INFO' column copied into the vep file
        :param vep_in: location of the vep file that needs extra information added to it.
        :param out_file: The location and name of where the output is supposed to be written to.
        """
        self.vcf_in = vcf_in
        self.vep_in = vep_in
        self.out_file = out_file

    def protocol(self):
        """
        Function that runs the protocol for this class.
        1. Create 2 objects of VcfAnalyser one for each input file.
        2. Add the columns to the VcfAnalyser object of the vep file.
        3. Add to the vep header to make it contain all needed information.
        4. Write the output to the requested output file.
        """
        # 1.
        self.vcf_tsv = vcfa.VcfAnalyser(self.vcf_in, "ID")
        self.vep_tsv = vcfa.VcfAnalyser(self.vep_in)
        # 2.
        self.add_information()
        # 3.
        self.add_to_header()
        # 4.
        self.vep_tsv.vep_to_file(self.out_file, "Extra")

    def add_information(self):
        """
        Function that adds and corrects information from the vep file. For exact detail see class description
        :return:
        """
        self.vep_tsv.column_header.append("INFO")
        for row_key in self.vep_tsv.keys():
            vep_ID = self.vep_tsv[row_key]["Uploaded_variation"]
            self.check_insertion_info(self.vep_tsv[row_key])
            self.check_start_codon(self.vep_tsv[row_key])
            self.check_transcript_amplification(self.vep_tsv[row_key])
            self.check_transcript_ablations(self.vep_tsv[row_key])
            try:
                vcf_info = self.vcf_tsv[vep_ID]["INFO"]
                self.vep_tsv[row_key]["INFO"] = vcf_info
            except KeyError:
                #if an ID from the vep file is not present in the vcf. This happens only if vep could not read a certain
                #ID from the vcf file or the vcf file is not the same as was used to make the vep file.
                print("WARNING: ID: {} was not found in the vcf file. This can mean the ID was not recocnized by vep"+\
                      "or in case of many of these warnings the wrong vcf file was given.")

    def check_insertion_info(self, row_dict):
        """
        Function that checks if an insertion would cause a frame shift or inframe insertion.
        :param row_dict: Dictionary containing a row from the vep file.
        """
        if row_dict["Allele"] == "insertion" and "coding_sequence_variant" in row_dict["Consequence"]:
            insertion = self.vcf_tsv[row_dict["Uploaded_variation"]]["ALT"]
            #the lenght of the insertion -1 because the sequence is annotated with a 'n' at the start
            if insertion != "<INS>" and (len(insertion) - 1) % 3 != 0:
                row_dict["Consequence"] += ",frameshift_variant"
                # replace the consequence with high.
                row_dict["Extra"] = row_dict["Extra"].replace("MODIFIER", "HIGH").replace("LOW", "HIGH").replace(
                    "MODERATE", "HIGH")
            elif insertion != "<INS>" and (len(insertion) - 1) % 3 == 0:
                row_dict["Consequence"] += ",inframe_insertion"
                # replace the consequence with moderate --> this is what the vep documentation states.
                row_dict["Extra"] = row_dict["Extra"].replace("MODIFIER", "MODERATE").replace("LOW", "MODERATE")

    def check_start_codon(self, row_dict):
        """
        Function that checks if a start_lost should have been annotated incase of a deletion, or shoud not have been
        annotated incase of a duplication.
        :param row_dict: Dictionary containing a row from the vep file.
        """
        if row_dict["Allele"] == "deletion" and "5_prime_UTR_variant" in row_dict["Consequence"] and \
                ("coding_sequence_variant" or "frameshift_variant") in row_dict["Consequence"] and "start_lost" not in \
                row_dict["Consequence"]:
            row_dict["Consequence"] += ",start_lost"
            row_dict["Extra"] = row_dict["Extra"].replace("MODIFIER", "HIGH").replace("LOW", "HIGH").replace(
                "MODERATE", "HIGH")
        elif row_dict["Allele"] == "duplication" and "start_lost" in row_dict["Consequence"]:
            row_dict["Consequence"] = row_dict["Consequence"].replace("start_lost", "").replace(",start_retained_variant", "")
            
    def check_transcript_amplification(self, row_dict):
        """
        Function that checks if a certain duplication should be a transcript amplification. This is the case if
        the full gene is overlapped by the duplication.
        :param row_dict: Dictionary containing a row from the vep file.
        """
        if row_dict["Allele"] == "duplication" and "OverlapPC=100" in row_dict["Extra"] and \
        "transcript_amplification" not in row_dict["Consequence"]:
            row_dict["Consequence"] = "transcript_amplification"
            
    def check_transcript_ablations(self, row_dict):
        """
        Function that adds missing transcript ablations. This are deletions that remove a complete gene. These are 
        not consequenly annotated.
	:param row_dict: Dictionary containing a row from the vep file.
        """
        if row_dict["Allele"] == "deletion" and "OverlapPC=100" in row_dict["Extra"] and \
        "transcript_ablation" not in row_dict["Consequence"]:
            row_dict["Consequence"] = "transcript_ablation"
    
            

    def add_to_header(self):
        """
        Function that adds info from the vcf_header to the vep_header that contain information about the INFO column.
        """
        for line in self.vcf_tsv.header:
            if line.startswith("##INFO"):
                self.vep_tsv.header.append(line)


if __name__ == "__main__":
    # add check to see if correct file types
    args = get_arguments()
    AddInformationVep(args.vcf, args.vep, args.output).protocol()
