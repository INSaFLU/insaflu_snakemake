# INSaFLU Snakemake

## Install some system dependencies

   ```bash
   sudo apt install git build-essential zlib1g-dev
   ```

## Create Environment

### Install Mamba Forge:

   ```bash
   curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
   ```

   ```bash
   bash Mambaforge-Linux-x86_64.sh
   ```

### Activate Conda:

   ```bash
   conda activate base
   ```

### Create Environment:

   ```bash
   mamba env create --name insaflu --file config/insaflu.yaml
   ```

## Activate Environment
   ```bash 
   conda activate insaflu
   ```

## Setup the folder to run

### Change directory to user, create reads directory and move your reads to this directory
   ```bash
   cd user
   mkdir reads
   cp path/to/your/files/your_files.fastq.gz path/user/reads/
   cd ..
   ```
### Edit the csv in the folder user called project_metadata.csv
   - This csv has to have the following format:
      ```csv 
      sample_name,fastq1,fastq2,tech
      
      sample_01,sample_01_R1_001,sample_02_R2_002,illumina
      sample_02,samle_02,,ont
      ```
   - All samples much be in fastq.gz format (do not append .fastq.gz in the fastq1 and fastq2 fields)
### Go to user/project.yaml and fill the following fields
 - only_samples (Run the sample only pipeline or run like a project)
 - project_name (your project name)
 - fasta_reference (mandatory for project)
 - gb_reference (mandatory for project)
 - primers (mandatory for iVar project)
 - consensus_assembler (mandatory for project -> between "snippy" and "iVar")
 - abricate (soon will change to classification)

### Change the parameters in the run in user/parameters.yaml
   - Here is possible to change the parameters in some softwares, change at your own responsability. A copy of the default parameters are in config/default_parameters.yaml
   - Each project will recive a copy of this file inside projects/project_name_folder/
### The references are in user/references folder
To add new ones just copy your fasta and genbank file.
### You must keep the fasta identifier equal to the GenBank Locus
   ">ON563414" == "LOCUS ON563414"


### To change the names of the project.yaml, samples.csv and parameters.yaml it is in config/constants.yaml
It is possible to change the following entries:
   - sample_yaml
   - sample_csv
   - software_parameters
You should give the absolute path or the relative path starting in the results folder.

## Running the workflow

### Run your project or analyze your samples with the following command:
   ```bash
    snakemake -c {threads} --use-conda
   ```
