#!/bin/bash 
set -x
mkdir $CONDA_PREFIX/db/insaflu/
cp ../workflow/db/db_influenza_typing_v10.fasta $CONDA_PREFIX/db/insaflu/sequences
conda activate $CONDA_PREFIX
abricate --setupdb
