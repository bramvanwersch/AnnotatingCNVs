# AnnotatingCNVs
PLANNOTATOR is a pipeline for annotating plant CNVs using VEP(https://www.ensembl.org/info/docs/tools/vep/index.html) and visualising the resulting table in a dash(https://dash.plot.ly/) application aswel as an easy way of running Ontologizer(http://ontologizer.de/) for information on gene ontology terms.
### Installation 

##### Innitial installation
1. Create a folder for the plannotator installation:
```shell
mkdir plannotator
cd plannotator
```
2. Download the plannotator github:
```shell
git clone https://github.com/bramvanwersch/AnnotatingCNVs
```
3. Create a folder called ensembl-vep and go into it for the vep installation:
```shell
mkdir ensembl-vep
cd ensembl-vep
```
4. Run the install.py file in the AnnotatingCNVs/python_scripts folder. This will install VEP Ontologizer and nextflow and create a bin file for all the relevant scripts and programs. Vep is going to prompt for some input so read what is says and answer accordingly. Note make sure you are in the ensembl-vep folder or pearl packages get installed in the wrong place and VEP wil not work:
```shell
python ../AnnotatingCNVs/python_scripts/install.py
```
5. Add the following lines to your path or profile:
```shell
export PATH=$PATH:~/route/to/plannotator/folder/bin
export Ontologizer=~/route/to/plannotator/folder/Ontologizer.jar
```
6. To install the conda package(this assumes that conda is installed if this is not the case it can be downloaded here; https://docs.anaconda.com/anaconda/install/) that contains all python scripts use the envirionment.yaml file that was just downloaded from github use the commands below to get all modules installed in an evironment named plannotater.
```shell
conda config --add channels conda-forge
cd AnnotatingCNVs
conda env create -f environment.yaml
```
7. Activating the conda environment. This is something you will have to do every time you restart the terminal:
```shell
conda activate plannotator
```
##### Download a VEP cache
You will need to download some additional files to be able to run VEP:

1. Download a cache for vep if you want to. This depends on the species and if you have your own refernce genome or not. ftp://ftp.ensemblgenomes.org/pub/plants/current/variation/vep/ this is the place where you can find the latest releases. Make sure you are in the ensembl-vep folder. The code below is an example for Arabidopsis thaliana. 
```shell
cd caches
curl -O ftp://ftp.ensemblgenomes.org/pub/plants/current/variation/vep/arabidopsis_thaliana_vep_44_TAIR10.tar.gz
tar -xzf arabidopsis_thaliana_vep_44_TAIR10.tar.gz
../convert_cache.pl -d route/to/caches/folder --species arabidopsis_thaliana --version 44
```

##### Download Ontologizer files
For ontologizer there are a couple of files that need to be downloaded if you want to use Ontologizer.
1. Download the Obo file. This file can simply be downloaded from the geneontology website using:
```shell
wget http://purl.obolibrary.org/obo/go/go-basic.obo
```  
2. Downloading an association file(optional). This is a file that holds the association between the gene names and ontologizer terms, it is different for each organism and a list of available association files can be found on http://current.geneontology.org/products/pages/downloads.html.

3. Downloading the population file(optional). Is a file that lists all the genes of an organism seperated by newlines. Finding these files is different for each organism. Make sure that there are gene identifiers in these files and not transcript identifiers. If the file is filled with transcript identifiers you can use the bash command below to remove everything after the . and return all unique genes into a file named geneIdentifiers.txt
```shell
grep -oP ".+(?=\.)" nameOfFile | sort | uniq > geneIdentifiers.txt
``` 
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

Example of a command to run the annotation of your vcf file:
```shell
annotate_cnvs.nf --vcf location/of/vcf/file --cache_dir directory/of/downloaded/cache --cache_version cache/verion --species name/of/species --output_dir name/of/output/directory --obo obo/file --association association/file/of/species --population all/known/genes/of/species
```
Example of a command for running the dash application:
```shell
visualise_vep.py name/of/output/directory/of/annotation/vcf
```
