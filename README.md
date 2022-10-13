# insaflu_snakemake on linux

## Create environment 

### install mamba forge:
    
   $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
    
   $ bash Mambaforge-Linux-x86_64.sh

### activate conda:
    
   $ conda activate base

### create environment:
    
   $ mamba env create --name insaflu --file config/insaflu.yaml

## Activate Environment 
   $ conda activate insaflu

## Deactivate Enviroment
   $ conda deactivate insaflu

## Create new folder named 'user_data'
   $ mkdir user_data
   $ cp path/to/your/files/your_files path/to/user/data
## Create a csv in the folder config_user called sample_info.csv
 This csv has to have the following format: sample_name,fastq1,fastq2,tech
## Go to config_user.yaml and fill the fields
 - only_samples
 - project_name
 - fasta_reference
 - gb_reference
## Run your project or analyse your samples with the following command:
   $snakemake -c {threads} --use-conda