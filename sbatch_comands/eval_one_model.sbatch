#!/bin/bash
#SBATCH --job-name=scConcept_eval
#SBATCH --partition=dc-gpu-devel
#SBATCH --account=jinm16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --time=01:15:00
#SBATCH --output=/p/scratch/cjinm16/dipippo1/logs/%x_%j.out
#SBATCH --error=/p/scratch/cjinm16/dipippo1/logs/%x_%j.err

set -euo pipefail
#CHECKPOINT=$SCRATCH/dipippo1/scConcept/model_checkpoints/split_mouse_2/session_2/small_weighted__hgld4ax1/epochs/last.ckpt
#CHECKPOINT=$SCRATCH/dipippo1/scConcept/model_checkpoints/split_mouse/session_2/big_weighted__uv28i4o3/epochs/last.ckpt
CHECKPOINT=$SCRATCH/dipippo1/scConcept/model_checkpoints/finetune/3de8fprx/epochs/last.ckpt
echo "============================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Checkpoint: $CHECKPOINT"
echo "Start time: $(date)"
echo "============================================"

cd $HOME/projects
source activate_env.sh
cd scConcept-1/src

export CUDA_VISIBLE_DEVICES=0

#srun python -m concept.evaluate \
#  --checkpoint "$CHECKPOINT" \
#  --datasets isd zeng zhuang \
#  --umap \
#  --clisi cell_type \
#  --ilisi dataset \
#  --subsample 100000 


#srun python -m concept.evaluate -c "$CHECKPOINT" -d isd zeng zhuang --umap --pca 32
srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets isd zeng zhuang -s 1000000000000 --umap --clisi cell_type --ilisi dataset
#python -m concept.evaluate -c "$CHECKPOINT" -d atlas --atlas-panel all_all -s 50000 --umap
srun python -m concept.evaluate --checkpoint "$CHECKPOINT"  -d atlas --atlas-panel all --umap
#srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets atlas_train -s 10000000000000
#srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets atlas -s 1000000000000
#srun python -m concept.evaluate --checkpoint $CHECKPOINT --datasets isd zeng zhuang -s 100000000000 --umap


