#!/bin/bash
#SBATCH --job-name=scConcept_train
#SBATCH --partition=dc-gpu
#SBATCH --account=jinm16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:4
#SBATCH --mem=200G
#SBATCH --time=23:30:00
#SBATCH --output=/p/scratch/cjinm16/dipippo1/logs/%x_%j.out
#SBATCH --error=/p/scratch/cjinm16/dipippo1/logs/%x_%j.err

# Print job info
echo "============================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "GPUs: $SLURM_GPUS_ON_NODE"
echo "Start time: $(date)"
echo "============================================"

# Navigate to project
cd $HOME/projects

# Activate your environment
source activate_env.sh

# Navigate to the source directory
cd scConcept-1
export CUDA_VISIBLE_DEVICES=0,1,2,3
srun bash -c 'echo "task=$SLURM_PROCID gpu=$CUDA_VISIBLE_DEVICES"'

# Run training
echo "Starting training..."

# ============================================
# RUN CONFIGURATIONS (Hydra overrides)
# SCOMMENTA SOLO UNA srun ALLA VOLTA
# ============================================

echo "Starting training..."

# RUN 1️⃣ — BIG MODEL + UNIFORM + BATCH 256
#srun python src/concept/train.py split=split_mouse_2 wandb.run_name=big_uniform_batch256 datamodule.probabilistic_panel_sampling=False datamodule.dataloader.train.batch_size=256 model.dim_model=512 model.dim_hid=1024 model.nlayers=8

# RUN 2️⃣ — BIG MODEL + WEIGHTED + BATCH 256
#srun python src/concept/train.py split=split_mouse_2 wandb.run_name=big_weighted_batch256 datamodule.probabilistic_panel_sampling=True datamodule.dataloader.train.batch_size=256 model.dim_model=512 model.dim_hid=1024 model.nlayers=8

# RUN 3️⃣ — SMALL MODEL + UNIFORM + BATCH 256
#srun python src/concept/train.py split=split_mouse_2 wandb.run_name=small_uniform_batch256 datamodule.probabilistic_panel_sampling=False datamodule.dataloader.train.batch_size=256 model.dim_model=128 model.dim_hid=256 model.nlayers=6

# RUN 4️⃣ — SMALL MODEL + WEIGHTED + BATCH 256
srun python src/concept/train.py split=split_mouse_2 wandb.run_name=small_weighted_batch256 datamodule.probabilistic_panel_sampling=True datamodule.dataloader.train.batch_size=256 model.dim_model=128 model.dim_hid=256 model.nlayers=6

echo "============================================"
echo "End time: $(date)"
echo "============================================"

