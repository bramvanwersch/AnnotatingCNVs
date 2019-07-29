#!/home/wersc004/Data/Software/nextflow/nextflow

params.obo = ""
params.association = ""
params.population = ""
params.study = ""
params.output_dir

obo = file(params.obo)
association = file(params.association)
population = file(params.population)
study = file(params.study)
output = file(params.output_dir)

out_name =  study.baseName

process filter_population_genes{
	input:
	file population
	
	output:
	file "all_population_genes.txt" into all_population_genes
	
	"""
	grep -oP ".+(?=\\.)" ${population} | sort | uniq > all_population_genes.txt
	"""
}

process run_ontologizer{
	publishDir "${output}"	
	
	input:
	file obo
	file association
	file pop from all_population_genes
	file study
	
	output:
	file 'table-*' into significant_ontologizer_results
	file 'view-*' into ontologizer_image_file
	
	"""
	java -jar ~/Data/Software/ontologizer/Ontologizer.jar -g ${obo} \
	-a ${association} -p ${pop} -s ${study}\
	-m Benjamini-Hochberg -c Parent-Child-Union -d
	"""
}

process create_ontology_picture{
	publishDir "${output}", mode: "move"
	
	input:
	file image from ontologizer_image_file
	
	output:
	file '*.png' into go_graph_picture
	
	"""
	dot -Tpng ${image} -o ${out_name}.png
	"""
}
