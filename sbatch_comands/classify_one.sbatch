#!/bin/bash
#SBATCH --job-name=scConcept_eval
#SBATCH --partition=dc-gpu-devel
#SBATCH --account=jinm16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --time=01:45:00
#SBATCH --output=/p/scratch/cjinm16/dipippo1/logs/%x_%j.out
#SBATCH --error=/p/scratch/cjinm16/dipippo1/logs/%x_%j.err

set -euo pipefail

# ============================================
# CONFIGURATION - MODIFY THIS
# ============================================
#CHECKPOINT=$SCRATCH/dipippo1/scConcept/model_checkpoints/split_mouse_2/session_2/big_weighted__6iz4xrqm/epochs/last.ckpt

# ============================================
# SETUP
# ============================================
echo "============================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Checkpoint: $CHECKPOINT"
echo "Start time: $(date)"
echo "============================================"

cd $HOME/projects
source activate_env.sh
cd scConcept-1/src
export CUDA_VISIBLE_DEVICES=0

# ============================================
# EXTRACT MODEL NAME FROM CHECKPOINT
# ============================================
# Extract the model folder name (e.g., big_weighted__6iz4xrqm)
CKPT_DIR=$(dirname $(dirname $CHECKPOINT))
MODEL_FOLDER=$(basename $CKPT_DIR)
STEP_NAME=$(basename $CHECKPOINT .ckpt | sed 's/=//g')
MODEL_NAME="${MODEL_FOLDER}_${STEP_NAME}"

EMB_BASE="/p/home/jusers/dipippo1/jureca/projects/scConcept-1/src/concept/concept_embeddings/models/${MODEL_NAME}"

echo "Model name: $MODEL_NAME"
echo "Embeddings base: $EMB_BASE"

# ============================================
# CONFIGURATION
# ============================================
REMOVE_OLD_EMBEDDINGS=false  # Set to true to remove old embeddings

# ============================================
# REMOVE OLD EMBEDDINGS (if enabled)
# ============================================
if [ "$REMOVE_OLD_EMBEDDINGS" = true ]; then
    echo ""
    echo "============================================"
    echo "Removing old embeddings..."
    echo "============================================"
    for DATASET in atlas_train atlas zeng zhuang isd; do
        EMB_PATH="${EMB_BASE}/${DATASET}"
        if [ -d "$EMB_PATH" ]; then
            echo "Removing: $EMB_PATH"
            rm -rf "$EMB_PATH"
        fi
    done
else
    echo ""
    echo "============================================"
    echo "Skipping embedding removal (REMOVE_OLD_EMBEDDINGS=false)"
    echo "============================================"
fi

# ============================================
# GENERATE EMBEDDINGS
# ============================================
echo ""
echo "============================================"
echo "Generating embeddings..."
echo "============================================"

srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets atlas_train -s 1000000000000 || true
#srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets atlas -s 1000000000000 || true
#srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets zeng -s 1000000000000 || true
#srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets zhuang -s 1000000000000 || true
#srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets isd -s 1000000000000 || true

# ============================================
# CLASSIFICATION
# ============================================
echo ""
echo "============================================"
echo "Running classification..."
echo "============================================"

cd concept

# Atlas (pre-split)
srun python classify.py --dataset atlas --model $MODEL_NAME --label_key class || true

# Zeng, Zhuang, ISD (auto-split)
srun python classify.py u--dataset zeng --model $MODEL_NAME --label_key class || true
srun python classify.py --dataset zhuang --model $MODEL_NAME --label_key class || true
#srun python -m concept.classify --dataset isd --model $MODEL_NAME --label_key class || true

# Zeng, Zhuang, ISD - cell_type_mmc_raw
srun python classify.py --dataset zeng --model $MODEL_NAME --label_key cell_type_mmc_raw || true
srun python classify.py --dataset zhuang --model $MODEL_NAME --label_key cell_type_mmc_raw || true
srun python classify.py --dataset isd --model $MODEL_NAME --label_key cell_type_mmc_raw || true

echo ""
echo "============================================"
echo "Done! End time: $(date)"
echo "============================================"
