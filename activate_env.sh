#!/bin/bash

module purge
module load Stages/2025
module load GCC
module load Python
module load OpenMPI
module load matplotlib scikit-image scikit-learn JupyterLab git

VENV="$HOME/projects/jupyter/kernels/dipippo1_kernel"
source "$VENV/bin/activate"

echo "Environment activated:"
which python
python --version
