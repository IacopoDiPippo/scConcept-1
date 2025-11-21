import os

# --- W&B ---
os.environ["WANDB_MODE"] = "offline"
os.environ["WANDB_SILENT"] = "true"
os.environ["WANDB_CACHE_DIR"] = "/p/scratch/cjinm16/dipippo1/wandb/cache"
os.environ["WANDB_MEDIA_DIR"] = "/p/scratch/cjinm16/dipippo1/wandb/media"

os.environ["MASTER_ADDR"] = "127.0.0.1"
os.environ["MASTER_PORT"] = "29501"
os.environ["NCCL_IB_DISABLE"] = "1"
os.environ["NCCL_SOCKET_IFNAME"] = "eth0"
os.environ["GLOO_SOCKET_IFNAME"] = "eth0"

import sys
import shutil
import filecmp
import hydra
import h5py
import torch
import datetime
import lightning as L
import pandas as pd
from lightning.pytorch.callbacks import LearningRateMonitor, ModelCheckpoint
from lightning.pytorch.loggers import WandbLogger
from omegaconf import DictConfig, OmegaConf
from lightning.pytorch.utilities import rank_zero_only

from lamin_dataloader.dataset import GeneIdTokenizer
from concept.data.datamodules import AnnDataModule
from concept.data.utils import add_count_nnz
from concept.model import BiEncoderContrastiveModel
import wandb
from lightning.pytorch.strategies import DDPStrategy, ParallelStrategy, SingleDeviceStrategy
from pathlib import Path
# load hydra config without decorative:
from hydra import compose, initialize
# with initialize(version_base=None, config_path="conf", job_name="test_app"):
#     cfg = compose(config_name="config", overrides=overrides)
#     print(OmegaConf.to_yaml(cfg))



def train() -> None:
    
    bash_cfg = OmegaConf.from_cli()
    resume_from_checkpoint = bash_cfg.pop("resume_from_checkpoint", False)

    def is_wandb_offline():
        return os.environ.get("WANDB_MODE", "").lower() == "offline"

    if resume_from_checkpoint:

        checkpoint_path = bash_cfg.pop("checkpoint")

        print(f"Resuming from checkpoint: {checkpoint_path}")

        # -------------------------------------------
        # CASE 1: WANDB OFFLINE → Load config from checkpoint
        # -------------------------------------------
        if is_wandb_offline():
            print("W&B is OFFLINE → loading config from checkpoint metadata")
            ckpt = torch.load(checkpoint_path, map_location="cpu")

            # The config is stored inside Lightning's checkpoint
            if "hyper_parameters" in ckpt:
                cfg = OmegaConf.create(ckpt["hyper_parameters"]["config"])
            else:
                raise RuntimeError(
                    "Checkpoint does not contain embedded config. "
                    "You must save config into checkpoints by passing it as hyperparams."
                )

            # merge additional CLI overrides
            cfg = OmegaConf.merge(cfg, bash_cfg)
            print(OmegaConf.to_yaml(cfg))
        
        # -------------------------------------------
        # CASE 2: WANDB ONLINE → load config from W&B API
        # -------------------------------------------
        else:
            print("W&B ONLINE → loading config from wandb cloud")

            run_id = bash_cfg.pop("run_id")
            wandb.login()
            api = wandb.Api()
            run = api.run(f'{bash_cfg.wandb.entity}/{bash_cfg.wandb.project}/{run_id}')
            
            cfg = DictConfig(run.config)
            cfg = OmegaConf.merge(cfg, bash_cfg)
            print(OmegaConf.to_yaml(cfg))

    else:
        print("Starting new training ...")
        print('overrides:', sys.argv[1:])
        with initialize(version_base=None, config_path="./conf"):
            cfg = compose(config_name="config", overrides=sys.argv[1:])

    
    dataset_path = cfg.PATH.ADATA_PATH
    if cfg.PATH.LOCAL_DIR is not None:
        local_dir = cfg.PATH.LOCAL_DIR
        assert len(set(local_dir.split("/"))) > 2, f"local_dir should not be root directory: {local_dir}"
        print(f'Copying training anndata files to directory...')
        if not os.path.exists(local_dir):
            os.makedirs(local_dir, exist_ok=True)
        copy_count = 0
        for _, filenames in cfg.PATH.SPLIT.items():
            if filenames is None:
                continue
            for file in filenames:
                src_file = os.path.join(dataset_path, file)
                dst_file = os.path.join(local_dir, file)
                if not os.path.exists(dst_file) or not filecmp.cmp(src_file, dst_file):
                # if not os.path.exists(dst_file):
                    shutil.copy(src_file, dst_file)
                    copy_count += 1
        print(f'{copy_count} files copied successfully!')
        dataset_path = local_dir
        
        # for _, filenames in cfg.PATH.SPLIT.items():
        #     if filenames is None:
        #         continue
        #     for file in filenames:
        #         src_file = os.path.join(local_dir, file)
        #         with h5py.File(src_file) as f:
        #             if "var/_count_nnz" in f:
        #                 print(f'Count nnz already exists in {file}...')
        #                 continue
        #         print(f'Adding count nnz for {src_file}...')
                # add_count_nnz(src_file)
    
    
    if 'val' in cfg.datamodule.dataset and cfg.datamodule.dataset.val is not None:
        val_loader_names = sorted(list(cfg.datamodule.dataset.val.keys()))
    else:
        val_loader_names = []
        
    gene_mapping = pd.read_pickle(cfg.PATH.gene_mapping_path).to_dict()
    

    datamodule_args = {    
        'dataset_path': dataset_path,
        'split': cfg.PATH.SPLIT,
        'panels_path': cfg.PATH.PANELS_PATH,
        'columns': cfg.datamodule.columns,
        'precomp_embs_key': cfg.datamodule.precomp_embs_key,
        'normalization': cfg.datamodule.normalization,
        'gene_sampling_strategy': cfg.datamodule.gene_sampling_strategy,
        'model_speed_sanity_check': cfg.datamodule.model_speed_sanity_check,
        # make sure to pass a copy to vaoid being modified before uploading to wandb:
        'dataset_kwargs': {**OmegaConf.to_container(cfg.datamodule.dataset, resolve=True, throw_on_missing=True)}, 
        'dataloader_kwargs': {**OmegaConf.to_container(cfg.datamodule.dataloader, resolve=True, throw_on_missing=True)},
        'val_loader_names': val_loader_names,
        'tokenizer': GeneIdTokenizer(gene_mapping)
    }
    datamodule = AnnDataModule(**datamodule_args)

    if cfg.wandb.enabled:
        logger = WandbLogger(name=cfg.wandb.run_name, entity=cfg.wandb.entity, project=cfg.wandb.project, save_dir=cfg.PATH.SCRATCH_ROOT, log_model=False, mode="offline",)
    
    CHECKPOINT_PATH = "dummy"
    if rank_zero_only.rank == 0:
        CHECKPOINT_PATH = os.path.join(cfg.PATH.CHECKPOINT_ROOT, logger.experiment.id if cfg.wandb.enabled else 'dummy')
        if cfg.wandb.enabled:
            logger.experiment.config.update(OmegaConf.to_container(cfg, resolve=True, throw_on_missing=True))
        print(OmegaConf.to_yaml(cfg))

    LIGHTNING_LOG_DIR = os.path.join(cfg.PATH.SCRATCH_ROOT, "lightning_logs")
    os.makedirs(LIGHTNING_LOG_DIR, exist_ok=True)


    trainer_kwargs = {
        'max_steps': cfg.model.training.max_steps,
        'accelerator': cfg.model.training.accelerator,
        'devices': cfg.model.training.devices,
        'num_nodes': int(os.environ['SLURM_JOB_NUM_NODES']) if 'SLURM_JOB_NUM_NODES' in os.environ else cfg.model.training.num_nodes,
        'logger': logger if cfg.wandb.enabled else None,
        'val_check_interval': cfg.model.training.val_check_interval,
        'check_val_every_n_epoch': cfg.model.training.check_val_every_n_epoch,
        'limit_train_batches': cfg.model.training.limit_train_batches,
        'callbacks': [
            LearningRateMonitor(logging_interval='step'),
            ModelCheckpoint(dirpath=CHECKPOINT_PATH, filename='min_train_loss', monitor='train/loss', mode='min',
                            every_n_epochs=1, save_top_k=1),
            ModelCheckpoint(dirpath=CHECKPOINT_PATH, filename='min_val_loss', monitor='val/loss', mode='min',
                            every_n_epochs=1, save_top_k=1),
            ModelCheckpoint(dirpath=os.path.join(CHECKPOINT_PATH, 'epochs'), filename='{epoch}', every_n_epochs=1, save_on_train_epoch_end=True, save_top_k=-1, save_last='link'),
            ModelCheckpoint(dirpath=os.path.join(CHECKPOINT_PATH, 'steps'), filename='{step}', every_n_train_steps=10000, monitor='train/loss', save_top_k=-1), # save a checkpoint every 10K steps
        ],
    }
    trainer_kwargs["default_root_dir"] = LIGHTNING_LOG_DIR
    # Automatically pick the right strategy
    if int(cfg.model.training.devices) > 1 or int(cfg.model.training.num_nodes) > 1:
        strategy = DDPStrategy(find_unused_parameters=True)
    else:
        strategy = SingleDeviceStrategy(device=torch.device("cuda:0"))

    trainer = L.Trainer(**trainer_kwargs, 
                        strategy=strategy,
                        precision='bf16-mixed', 
                        use_distributed_sampler=False,
                        accumulate_grad_batches=cfg.model.training.accumulate_grad_batches,
                        )



    model_args = {
        'config': cfg.model,
        'pad_token_id': gene_mapping['<pad>'],
        'cls_token_id': gene_mapping['<cls>'],
        'vocab_size': len(gene_mapping),
        'world_size': trainer.world_size, 
        'val_loader_names': val_loader_names, 
        'precomp_embs_key': cfg.datamodule.precomp_embs_key,
    }
    model = BiEncoderContrastiveModel(**model_args)

    if not resume_from_checkpoint and cfg.model.training.validate_before_training:
        trainer.validate(model=model, 
                        datamodule=datamodule,
                        )
    
    if resume_from_checkpoint:
        ckpt_path = os.path.join(cfg.PATH.CHECKPOINT_ROOT, run_id, checkpoint)
        model = BiEncoderContrastiveModel.load_from_checkpoint(ckpt_path, **model_args, strict=False)
    
    trainer.fit(model=model, 
                datamodule = datamodule)

    # if cfg.PATH.LOCAL_DIR is not None and 'persistent' not in cfg.PATH.LOCAL_DIR:
    #     print(f'Removing trianing local directory...')
    #     shutil.rmtree(local_dir)
    #     print(f'Local directory removed!')

if __name__ == "__main__":
    train()
