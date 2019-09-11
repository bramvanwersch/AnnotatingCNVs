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
    #system("git clone https://github.com/Ensembl/ensembl-vep.git {}".format(location))
    #system("git checkout release/96 {}".format(location))
    #system("perl {}/INSTALL.pl".format(location))
    #system("mkdir {}/caches".format(location))
    system("wget http://ontologizer.de/cmdline/Ontologizer.jar -P ../{}".format(location))
    system("cd ../{} && curl -s https://get.nextflow.io | bash".format(location))
    system("mkdir -p ../{}/bin".format(location))
    for file_loc in program_file_list:
        system("ln -sf ../{}/{} ../{}/bin".format(location,file_loc, location))

if __name__ == "__main__":
    loc = argv[1]
    install(loc)
