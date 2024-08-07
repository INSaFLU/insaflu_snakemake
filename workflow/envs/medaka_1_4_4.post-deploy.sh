conda activate $CONDA_PREFIX
conda update bcftools --channel conda-forge
conda install -c conda-forge gsl

set -x
cp ../workflow/software/replacement_scripts/run_check_consensus $CONDA_PREFIX/bin/run_check_consensus
cp ../workflow/software/replacement_scripts/medaka_consensus $CONDA_PREFIX/bin/medaka_consensus
cp ../workflow/software/replacement_scripts/ivar $CONDA_PREFIX/bin/ivar
env_name="ivar8745920361"

if conda env list | grep -q "$env_name"; then
    echo "Environment $env_name already exists."
else
    mamba env create --name "$env_name" --file ../workflow/envs/ivar.yaml
fi

