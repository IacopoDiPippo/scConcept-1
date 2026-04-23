#!/bin/bash
#SBATCH --job-name=scConcept_train_noval_noisd
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


# RUN no val no isd
#srun python src/concept/train.py split=split_mouse2_noval  wandb.run_name=small_weighted_noisd datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.check_val_every_n_epoch=1000000 datamodule=NoValBasic

#srun python src/concept/train.py split=split_mouse2_noval  wandb.run_name=small_weighted_noisd datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.check_val_every_n_epoch=1000000 datamodule=NoValBasic initialize_2.resume=True initialize_2.run_id=sj7tx2hc initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=False

#srun python src/concept/train.py split=split_mouse2_noval  wandb.run_name=small_weighted_noisd datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.check_val_every_n_epoch=1000000 datamodule=NoValBasic initialize_2.resume=True initialize_2.run_id=ykbrnhcf initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=False


#srun python src/concept/train.py split=split_spatial_noisd_noval  wandb.run_name=small_weighted_noisd_zeng_zhang_finetuning  datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.check_val_every_n_epoch=1000000 datamodule=DataModuleSpatialNOVAL  datamodule.columns=[] datamodule.dataset.train.transcriptomics.panel_max_drop_rate=0.1 datamodule.dataset.train.transcriptomics.panel_selection_mixed_prob=0 datamodule.finetune_panels=True  initialize_2.resume=True initialize_2.run_id=yga0ei93 initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=False 


#srun python src/concept/train.py split=split_mouse2_noval  wandb.run_name=small_weighted_noisd datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.check_val_every_n_epoch=1000000 datamodule=NoValBasic initialize_2.resume=True initialize_2.run_id=yga0ei93 initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=False


#srun python src/concept/train.py split=split_spatial_noisd_noval  wandb.run_name=small_weighted_noisd_zeng_zhang_finetuning  datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.check_val_every_n_epoch=1000000 datamodule=DataModuleSpatialNOVAL  datamodule.columns=[] datamodule.dataset.train.transcriptomics.panel_max_drop_rate=0.1 datamodule.dataset.train.transcriptomics.panel_selection_mixed_prob=0 datamodule.finetune_panels=True  initialize_2.resume=True initialize_2.run_id=bme2g689 initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=False 

srun python src/concept/train.py split=split_spatial_noisd_noval  wandb.run_name=small_weighted_noisd_zeng_zhang_finetuning  datamodule.probabilistic_panel_sampling=True  model.dim_model=128 model.dim_hid=256 model.nlayers=6 model.training.check_val_every_n_epoch=1000000 datamodule=DataModuleSpatialNOVAL  datamodule.columns=[] datamodule.dataset.train.transcriptomics.panel_max_drop_rate=0.1 datamodule.dataset.train.transcriptomics.panel_selection_mixed_prob=0 datamodule.finetune_panels=True  initialize_2.resume=True initialize_2.run_id=mhvog2y2 initialize_2.checkpoint=epochs/last.ckpt initialize_2.do_validate_first=False
echo "============================================"
echo "End time: $(date)"
echo "============================================"

                                                                                                                                                                                                                                                    65,0-1        Bot