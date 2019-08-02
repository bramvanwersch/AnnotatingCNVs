#!/home/wersc004/Data/Software/nextflow/nextflow


/*Autor: Bram van Wersch
*University: Wageningen university
*/

params.vcf = ""
params.output_dir = ""
params.cache_dir = ""
params.cache_version = ""
params.species = ""
params.gff = ""
params.fasta = ""
params.run_ontology = true
params.obo = ""
params.association = ""
params.population = ""

/*FILE SPECIFIED CHECK
*check if ontologizer files are specified if not make sure to raise a custom
*error instead of the one that will be raised by nextflow. Besided that 
*make sure the user enters false if the ontologizer is not supposed to run
*to prevent mistakes
*/

if ((params.obo == "" || params.association == "" || params.population == "") && params.run_ontology == true)
	exit 1, "Please specify an --obo, --association and --population file"
	
vcf = file(params.vcf)

if (params.gff != "" && params.fasta != ""){
	gff_in = file(params.gff)
	fasta = file(params.fasta)
}

if (params.run_ontology == true){
	obo = file(params.obo)
	association = file(params.association)
	population = file(params.population)
    ontologizer_file_list = [params.obo, params.association, params.population]
    ontologizer_name_list = ["obo_file.obo", "association_file.gaf", "population_file.txt"]
}
else if (params.run_ontology != false){
	error "Run_ontology argument can be true or false not ${params.run_ontology}"
}

output = file(params.output_dir)

/*
*print out a starting sequence so people can verify their input files.
*/

log.info "\nANNOTING PLANT CNVs"
log.info "===================================="
log.info "VCF: ${params.vcf}"
log.info "output directory: ${params.output_dir}"
log.info "cache directory: ${params.cache_dir}"
log.info "cache version: ${params.cache_version}"
log.info "species: ${params.species}"
log.info "gff: ${params.gff}"
log.info "fasta: ${params.fasta}"
log.info "Run ontology: ${params.run_ontology}"
log.info "obo: ${params.obo}"
log.info "association: ${params.association}"
log.info "population: ${params.population}"

/*
*FILE AND INPUT CHECK.
*Check if a fasta and gff file are given. If neither was given it is 
*assumed the user wants to use the cache mode of VEP. If both are given 
*it is assumed the user wants to use the custom mode of VEP and if only
*one of them is given it is assumed the user made a mistake and an eror 
*will be printed.
*/

if (params.gff != "" && params.fasta != "")
	mode = "custom"
else if (params.gff == "" && params.fasta == ""){
	if (params.cache_version == "" || params.species == "" || params.cache_dir == ""){
		exit 1, "Missing one of the variables --species, --cache_dir or --cache_verion"
		}
	else{
		mode = "cache"
		}
	}
else{
	if (params.fasta == "") 
		exit 1, "Missing fasta file, cannot find given file: ${fasta}"
	if (params.gff == "") 
		exit 1, "Missing gff file, cannot find given file: ${gff}"
	}

if (params.run_ontology == "yes"){
	if (params.obo == "") 
		exit 1, "Missing obo file, cannot find given file: ${obo}"
	if (params.association == "") 
		exit 1, "Missing association file, cannot find given file: ${association}"
	if (params.poulation == "") 
		exit 1, "Missing population file, cannot find given file: ${population}"
}

log.info "Running in '${mode}' mode"
log.info "\n"

/*
*Adding insertions caused by dispersed duplications.
*/
process add_dispersed_insertions{

	input:
	file vcf
	
	output:
	file 'added_vcf.vcf' into added_vcf
	
	"""
	add_dispersed_duplications.py \
	-vcf ${vcf} --output added_vcf.vcf
	"""
}

/*
*Bgzip the file and index it with tabix if the user wants to run vep with 
*a gff and fasta file.
*/
if (mode == "custom"){
	process create_bgzip_file{
		publishDir "${output}/temp_dir"
		
		input:
		file gff_in
		
		output:
		file 'gff_file.gff.gz' into gff_bgzip
		
		"""
		grep -v "#" ${gff_in} | sort -k1,1 -k4,4n -k5,5n -t\$'\\t' |\
		bgzip -c > gff_file.gff.gz
		"""
	}

	process create_tabbix_file{
		publishDir "${output}/temp_dir"
		
		input:
		file gff_bgzip
		
		output:
		file '*.tbi' into gff_tbi
		
		"""
		tabix -f -p gff ${gff_bgzip}
		"""
	}
}

/*
*Run vep and depending on the mode use 2 different commands. After running
*remove the temp dir directory that contains the bgzipped and tabbix file
*because they had to be located here to be able to find them
*/
process run_vep{

	input:
	file input_file from added_vcf
	val species from params.species
	val cache_version from params.cache_version
	val cache from params.cache_dir
	
	output:
	file 'vep_output.txt' into vep_result
	
	afterScript "rm -r ${output}/temp_dir/"
	
	script:
	if (mode == "cache")
		"""
		vep --cache --offline --numbers --force_overwrite --dir ${cache}\
		--species ${species} --cache_version ${cache_version}\
		 -i ${input_file} -o vep_output.txt
		"""
		
	else if (mode == "custom")
		"""
		/mnt/scratch/timme068/ensembl-vep/vep --numbers --force_overwrite \
		--gff ${output}/temp_dir/gff_file.gff.gz --fasta ${fasta}\
		 -i ${input_file} -o vep_output.txt
		"""
		
	else
		error "Invalid allignment mode: ${mode}. Choose from custom or cache(default)"
	
}

/*
*Add some information that is lacking to the vep file for a more complete
*picture
*/
process add_info_vep_file{
	publishDir "${output}"	

	input:
	file vep_input_file from vep_result
	file vcf_input_file from added_vcf
	
	output:
	file 'added_vep_output.txt' into added_vep_result

	"""
	vep_add_info.py \
	--vcf ${vcf_input_file} --vep ${vep_input_file} \
	--output added_vep_output.txt
	"""
}

/*
*If the user requested to run ontologizer and all parameters are present
*run the functions required for ontologizer
*/
if (params.run_ontology == true){
    /*
    *Get all genes that are in deletions and overlapping a coding region
    *as test set for the ontologizer run.
    */
	process get_deletion_gene_names{
		
		input:
		file input_file from added_vep_result
		
		output:
		file 'all_deletion_coding_cnvs.txt' into all_deletion_identifiers
		
		"""
		grep -vP "^#" ${input_file}\
		 | grep -P "deletion.+(coding_sequence_variant|frameshift_variant)"\
		 | cut -f 4 | sort |uniq > all_deletion_coding_cnvs.txt
		"""
	}
    /*
    *Filter the list of population genes on transcripts and make sure 
    *only gene names are actualy present
    */
	process filter_population_genes{

		input:
		file population
		
		output:
		file "all_population_genes.txt" into all_population_genes
		
		"""
		grep -oP ".+(?=\\.)" ${population} | sort | uniq > all_population_genes.txt
		"""
	}
    
    /*
    *Run ontologizer if the study file is not empty
    */
	process run_ontologizer{

		publishDir "${output}"
            
		input:
		file obo
		file association
		file pop from all_population_genes
		file study from all_deletion_identifiers
		
		output:
		file 'table-*' into significant_ontologizer_results
		file 'view-*' into ontologizer_image_file
        
        """
        java -jar $Ontologizer -g ${obo} \
        -a ${association} -p ${pop} -s ${study}\
        -m Benjamini-Hochberg -c Parent-Child-Union -d
        """
	}
    
    /*
    *Create symbolic links to the ontologizer files to make sure the 
    *dash script can find them later on. Also make sure to create the
    *directory if it does not already exist.
    */
    process create_symbolic_links{
        input:
        val loc from ontologizer_file_list
        val name from ontologizer_name_list

        beforeScript "mkdir -p ${output}"

        """
        ln -sf ${loc} ${output}/${name}
        """
    }

    /*
    *Create a picture using dot if the view file was created by ontologizer
    */
	process create_ontology_picture{

		publishDir "${output}", mode: "move"
		
		input:
		file image from ontologizer_image_file
		
		output:
		file "GO_graph.png" into go_graph_picture
		
		"""
		dot -Tpng ${image} -o GO_graph.png
		""" 
	}
}

log.info "\n##########VISUALISE##########"
log.info "Use the following command to run a dash app that visualised \
the vep and ontologizer results: visualise_vep.py ${output}"
log.info "\n"


