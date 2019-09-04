# AnnotatingCNVs
PLANNOTATOR is a pipeline for annotating plant CNVs using VEP(link to VEP) and visualising the resulting table in a dash application aswel as an easy way of running Ontologizer(link to ontologizer) for information on gene ontology terms.
### Installation 
1. Create a folder to download all the programs except for the conda package into:
```shell
mkdir some_name
cd some_name
``` 
2. Download VEP version 96 and make sure you download HTSLIB if you want to use a custom annotation genome.
```shell
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git pull
git checkout release/96
perl INSTALL.pl
```
3. Download a cache for vep if you want to. This depends on the species and if you have your own refernce genome or not. ftp://ftp.ensemblgenomes.org/pub/plants/current/variation/vep/ this is the place where you can find the latest release. Make sure you are in the ensembl_vep folder. This should have create a cache folder with a folder named after your species in it that contains the ensembl-cache. 
```shell
mkdir caches
curl -O ftp://ftp.ensemblgenomes.org/pub/plants/current/variation/vep/name_of_the_desired_cache
tar -xzf arabidopsis_thaliana_vep_43_TAIR10.tar.gz
convert_cache.pl -d route/to/VEP/installation/ ensembl-vep/caches --species your_species --version version_off_cache
```
4. Download Ontologizer. Make sure you are in the some_name folder.
```shell

```
### Getting the required files
Here is a short description for getting most of the files that you will need for running ontologizer. If these files are not available from the sources noted below then there is a chance you can find them somewhere on the internet but in my experience that chance is not to high. Alternatives are to create your own using available tools or use a webtool that has these files simply available. Study sets can be downloaded after running the final script to run your webtool with.

__Obo file__: wget http://purl.obolibrary.org/obo/go/go-basic.obo

__Association file__: is different for each organism and not all organisms have an association file. List of available files can be found on http://current.geneontology.org/products/pages/downloads.html

__Population file__: Is a file that lists all the genes of an organism seperated by newlines. Finding these files is different for each organism. Make sure that there are gene identifiers in these files and not transcript identifiers. If the file is filled with transcript identifiers you can use:
- grep -oP ".+(?=\.)" nameOfFile | sort | uniq > geneIdentifiers.txt

This removes the .number extension of all the transcripts and filters out all the double names and puts it into a file called geneIdentifiers.txt

### Usage
To start the pipeline run the annotate_cnvs.nf script. This is the main pipeline script, there are allot necesairy files to run the full pipeline. Here is a list of them:

Command | Required | Explanation | default 
--- | --- | --- |---
--vcf | yes | file that contains your CNVs in a vcf format.
--output_dir | yes | name of the directory the output will be created in. Be carefull will overwrite existing directories.
--cache_dir | if running VEP with cache | directory containing the caches of the organisms
--cache_version | if running VEP with cache | version of the downloaded cache see also VEP installation
--species | if running VEP with cache | name of the species. This is zo VEP can find the cache.
--gff | if running VEP with custom annotation | gff file for the organism of intrest.
--fasta | if running VEP with custom annotation | fasta file for the organism of intrest.
--run_ontology | no | Either true or false based depending if you want to run ontologizer | true
--obo | if run_ontology = true | File containing the ontology graph and the connections between them. See getting files section for more information
--association | if run_ontology = true | File that contains the associations of ontology terms with the genes of an organism. See getting files section for more information
--population  | if run_ontology = true | All genes (not transcripts) of the organism of interest. See getting files section for more information

Example of a command 
annotate_cnvs.nf --vcf location/of/vcf/file --cache_dir directory/of/downloaded/cache --cache_version cache/verion --species name/of/species --output_dir name/of/output/directory --obo obo/file --association association/file/of/species --population all/known/genes/of/species
