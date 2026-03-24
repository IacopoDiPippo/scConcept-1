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

cd $HOME/projects
source activate_env.sh
cd scConcept-1

export CUDA_VISIBLE_DEVICES=0,1,2,3
srun bash -c 'echo "task=$SLURM_PROCID gpu=$CUDA_VISIBLE_DEVICES"'

# ============================================
# SCOMMENTA SOLO UNA srun ALLA VOLTA
# bs=512, lr scaling: 4e-4 base, x3 per ogni step
# ============================================

echo "Starting training..."

# lr=4e-4
#srun python src/concept/train.py wandb.run_name=lr4e-4 split=split_mouse_2 datamodule.probabilistic_panel_sampling=True datamodule.dataloader.train.batch_size=512 model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.lr=4e-4 model.training.max_steps=250000 model.training.check_val_every_n_epoch=5
srun python src/concept/train.py wandb.run_name=little_output_space split=split_mouse_2   datamodule.probabilistic_panel_sampling=True model.projection_dim=32  datamodule.columns='[]'  datamodule.dataloader.train.transcriptomics.batch_size=512 model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.lr=12e-4  model.training.check_val_every_n_epoch=5
# lr=12e-4
#srun python src/concept/train.py wandb.run_name=little_output_space split=split_mouse_2   datamodule.probabilistic_panel_sampling=True model.projection_dim=32  datamodule.dataloader.train.transcriptomics.batch_size=512 model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.lr=12e-4  model.training.check_val_every_n_epoch=5

# lr=36e-4
#srun python src/concept/train.py wandb.run_name=lr36e-4 split=split_mouse_2 datamodule.probabilistic_panel_sampling=True datamodule.dataloader.train.batch_size=512 model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.lr=36e-4 model.training.max_steps=250000 model.training.check_val_every_n_epoch=5

# lr=108e-4
#srun python src/concept/train.py wandb.run_name=lr108e-4 split=split_mouse_2 datamodule.probabilistic_panel_sampling=True datamodule.dataloader.train.batch_size=512 model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.lr=108e-4 model.training.max_steps=250000 model.training.check_val_every_n_epoch=5

echo "Training finished at: $(date)"
