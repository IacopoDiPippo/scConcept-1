#!/bin/bash
#SBATCH --job-name=scConcept_train
#SBATCH --partition=dc-gpu-devel
#SBATCH --account=jinm16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:4
#SBATCH --mem=200G
#SBATCH --time=01:15:00
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
srun python src/concept/train.py wandb.run_name=session3spatial_small_weighted_finetune datamodule.columns='[]' split=split_spatial datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 initialize_2.resume=True initialize_2.run_id=kdp2po4o initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=True datamodule.dataset.train.transcriptomics.panel_max_drop_rate=0.1 datamodule.dataset.train.transcriptomics.panel_selection_mixed_prob=0 datamodule.finetune_panels=True  initialize_2.split=split_mouse_ISD datamodule=DataModuleSpatial

#srun python src/concept/train.py
#srun python -u  src/concept/train.py wandb.run_name=big_uniform_batch256 datamodule.probabilistic_panel_sampling=False datamodule.dataloader.train.batch_size=256 model.dim_model=512 model.dim_hid=1024 model.nlayers=8 split=split_mouse_2  model.training.limit_train_batches=300
#srun python src/concept/train.py wandb.run_name=little_output_space split=split_mouse_2   datamodule.probabilistic_panel_sampling=True model.projection_dim=32  datamodule.columns='[]'  datamodule.dataloader.train.transcriptomics.batch_size=512 model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.lr=12e-4  model.training.check_val_every_n_epoch=5
#srun python src/concept/train.py wandb.run_name=session3spatial_small_weighted datamodule.columns=[] split=split_only_spatial datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 initialize_2.resume=True initialize_2.run_id=kdp2po4o initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=True datamodule.dataset.train.transcriptomics.panel_max_drop_rate=0.1 datamodule.dataset.train.transcriptomics.panel_selection_mixed_prob=0 datamodule.finetune_panels=False  initialize_2.split=split_mouse
#srun python src/concept/train.py wandb.run_name=session3spatial_small_weighted datamodule.columns='[]'  split=split_spatial datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 initialize_2.resume=True initialize_2.run_id=kdp2po4o initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=True datamodule.dataset.train.transcriptomics.panel_max_drop_rate=0.1 datamodule.dataset.train.transcriptomics.panel_selection_mixed_prob=0 datamodule.finetune_panels=False  initialize_2.split=split_mouse

#srun python src/concept/train.py wandb.run_name=session3P_small_weighted split=split_mouse datamodule.probabilistic_panel_sampling=True datamodule.dataloader.train.batch_size=256 model.dim_model=128 model.dim_hid=256 model.nlayers=6 initialize_2.resume=True initialize_2.run_id=kdp2po4o initialize_2.split=split_mouse initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=True datamodule.dataset.train.panel_max_drop_rate=0.1 datamodule.dataset.train.panel_selection_mixed_prob=0 datamodule.finetune_panels=True datamodule.dataset.train.panel_size_max=1200
echo "============================================"
echo "End time: $(date)"
echo "============================================"
