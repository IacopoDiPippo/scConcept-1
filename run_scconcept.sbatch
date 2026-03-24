#!/bin/bash
#SBATCH -J scconcept
#SBATCH --output=slurm-outputs/slurm-%j.out
#SBATCH --error=slurm-outputs/slurm-%j.out
#SBATCH -p dc-gpu
#SBATCH --account=jinm16
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --mem=200G
#SBATCH -t 01:30:00

source "$HOME/.bash_profile"

cd $HOME/projects
source activate_env.sh

export WANDB__SERVICE_WAIT=300
export WANDB_TAGS="slurm-$SLURM_JOB_ID"

echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

cd $HOME/projects/scConcept-1
srun python src/concept/train.py


