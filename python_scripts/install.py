#!/usr/bin/env python3

"""
Author: Bram van Wersch
University: Wageningen university
Date: 13/06/2019
"""

from sys import argv
from os import system

def install(location):
    program_file_list = ["AnnotatingCNVs/python_scripts/add_dispersed_insertions.py",\
                    "AnnotatingCNVs/python_scripts/correct_vep.py",\
                    "AnnotatingCNVs/python_scripts/vcf_analyser.py",\
                    "AnnotatingCNVs/python_scripts/visualise_vep.py",\
                    "AnnotatingCNVs/nextflow_scripts/annotate_cnvs.nf",\
                    "AnnotatingCNVs/nextflow_scripts/get_go_terms.nf",
                    "Ontologizer.jar",\
                    "ensembl-vep/vep",\
                    "nextflow"]
    system("git clone https://github.com/Ensembl/ensembl-vep.git .")
    #system("git checkout release/96 {}".format(location))
    #system("perl {}/INSTALL.pl".format(location))
    #system("mkdir {}/caches".format(location))
    system("wget http://ontologizer.de/cmdline/Ontologizer.jar -P ..")
    system("cd ../ && curl -s https://get.nextflow.io | bash")
    system("mkdir -p ../bin")
    for file_loc in program_file_list:
        system("ln -sf ../{} ../bin".format(file_loc))

if __name__ == "__main__":
    loc = argv[1]
    install(loc)
