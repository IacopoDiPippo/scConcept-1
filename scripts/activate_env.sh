#!/bin/bash
# ============================================================
#  Environment setup script for scConcept-1 project
#  Loads modules, activates virtual environment, sets paths
# ============================================================

# --- 1. Load required modules ---
module purge
module load Stages/2025
module load GCC
module load Python
module load OpenMPI
module load matplotlib scikit-image scikit-learn JupyterLab git

# --- 2. Get script directory (absolute path) ---
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# --- 3. Define virtual environment path ---
VENV_PATH="${SCRIPT_DIR}/../.env/scConcept-1"

# --- 4. Activate the virtual environment ---
if [ -d "${VENV_PATH}" ]; then
    source "${VENV_PATH}/bin/activate"
else
    echo "❌ Virtual environment not found at: ${VENV_PATH}"
    echo "Please create it first using: python -m venv ${VENV_PATH}"
    exit 1
fi

# --- 5. Verify activation ---
if [[ -n "$VIRTUAL_ENV" ]]; then
    echo "✅ Virtual environment 'scConcept-1' activated successfully."
else
    echo "❌ Failed to activate virtual environment."
    exit 1
fi

# --- 6. Set Python path (automatically detect Python version) ---
PYV=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
export PYTHONPATH=${VIRTUAL_ENV}/lib/python${PYV}/site-packages:${PYTHONPATH}

# --- 7. Confirmation ---
echo "✅ Environment successfully loaded and configured!"
echo "   Virtual env: ${VIRTUAL_ENV}"
echo "   Python path: $(which python)"
