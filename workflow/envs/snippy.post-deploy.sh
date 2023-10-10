#!/bin/bash 
set -x
cp ../workflow/software/replacement_scripts/run_check_consensus $CONDA_PREFIX/bin/run_check_consensus
cp ../workflow/software/replacement_scripts/snippy $CONDA_PREFIX/bin/snippy
cp ../workflow/software/replacement_scripts/ivar $CONDA_PREFIX/bin/ivar
env_name="ivar8745920361"

if conda env list | grep -q "$env_name"; then
    echo "Environment $env_name already exists."
else
    mamba env create --name "$env_name" --file ../workflow/envs/ivar.yaml
fi
