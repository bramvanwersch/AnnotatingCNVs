
params.location = ""
params.dispersed_insertions = "AnnotatingCNVs/python_scripts/add_dispersed_insertions.py"
output = file(params.location)

program_file_list2 = [params.dispersed_insertions]
program_file_list = ["AnnotatingCNVs/python_scripts/add_dispersed_insertions.py",\
                    "AnnotatingCNVs/python_scripts/correct_vep.py",\
                    "AnnotatingCNVs/python_scripts/vcf_analyser.py",\
                    "AnnotatingCNVs/python_scripts/visualise_vep.py",\
                    "AnnotatingCNVs/nextflow_scripts/annotate_cnvs.nf",\
                    "AnnotatingCNVs/nextflow_scripts/get_go_terms.nf",
                    "Ontologizer.jar",\
                    "ensembl-vep/vep"]
    
process fill_bin_folder{
    input:
    val loc from program_file_list

    beforeScript "mkdir -p ${output}/bin"

    """
    ln -sf ${output}/bin/${loc} 
    """
    }
