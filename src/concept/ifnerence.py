import os
import torch
import lightning as L
import anndata as ad
import pandas as pd
from hydra import compose, initialize
from omegaconf import OmegaConf

from concept.model import BiEncoderContrastiveModel   # UPDATE THIS
from lamin_dataloader.dataset import GeneIdTokenizer
from concept.data.datamodules import AnnDataModule                # UPDATE THIS


def load_config(overrides=None):
    if overrides is None:
        overrides = []
    with initialize(version_base=None, config_path="conf"):
        cfg = compose(config_name="config", overrides=overrides)
    return cfg


def load_model(cfg, ckpt_path, gene_mapping):
    model_args = {
        "config": cfg.model,
        "pad_token_id": gene_mapping["<pad>"],
        "cls_token_id": gene_mapping["<cls>"],
        "vocab_size": len(gene_mapping),
        "world_size": 1,
        "val_loader_names": [],
        "precomp_embs_key": cfg.datamodule.precomp_embs_key,
    }

    print(f"🔍 Loading checkpoint:\n{ckpt_path}")
    model = BiEncoderContrastiveModel.load_from_checkpoint(
        ckpt_path, **model_args, strict=False
    )
    return model


def run_prediction(model, datamodule):
    trainer = L.Trainer(accelerator="gpu", devices=1, logger=False)
    print("\n🚀 Running prediction ...")

    preds = trainer.predict(model, datamodule=datamodule)

    print("🔄 Gathering embeddings ...")
    cls_emb = torch.cat([p["cls_cell_emb"] for p in preds], dim=0).cpu().numpy()
    mean_emb = torch.cat([p["mean_cell_emb"] for p in preds], dim=0).cpu().numpy()

    return cls_emb, mean_emb


def save_embeddings(adata_path, cls_emb, mean_emb):
    print(f"\n📂 Loading validation AnnData:\n{adata_path}")
    adata = ad.read_h5ad(adata_path)

    assert adata.n_obs == cls_emb.shape[0], (
        f"Cells mismatch: AnnData has {adata.n_obs}, embeddings have {cls_emb.shape[0]}"
    )

    adata.obsm["concept_cls_embedding"] = cls_emb
    adata.obsm["concept_mean_embedding"] = mean_emb

    out_path = adata_path.replace(".h5ad", "_concept_embeddings.h5ad")

    print(f"💾 Saving embeddings to:\n{out_path}")
    adata.write(out_path)

    print("🎉 Done!")


def main():

    # ------------------------------------------------
    # LOAD HYDRA CONFIG
    # ------------------------------------------------
    cfg = load_config()

    # mapping
    gene_mapping = pd.read_pickle(cfg.PATH.gene_mapping_path).to_dict()

    # ------------------------------------------------
    # BUILD DATAMODULE — *only val loader is used*
    # ------------------------------------------------
    datamodule_args = {
        "dataset_path": cfg.PATH.ADATA_PATH,
        "split": cfg.PATH.SPLIT,
        "panels_path": cfg.PATH.PANELS_PATH,
        "columns": cfg.datamodule.columns,
        "precomp_embs_key": cfg.datamodule.precomp_embs_key,
        "normalization": cfg.datamodule.normalization,
        "gene_sampling_strategy": cfg.datamodule.gene_sampling_strategy,
        "model_speed_sanity_check": False,
        "dataset_kwargs": OmegaConf.to_container(cfg.datamodule.dataset, resolve=True),
        "dataloader_kwargs": OmegaConf.to_container(cfg.datamodule.dataloader, resolve=True),
        "val_loader_names": list(cfg.datamodule.dataset.val.keys()),
        "tokenizer": GeneIdTokenizer(gene_mapping),
    }

    datamodule = AnnDataModule(**datamodule_args)

    # ------------------------------------------------
    # LOAD MODEL CHECKPOINT
    # ------------------------------------------------
    ckpt_path = os.path.join(
        cfg.PATH.CHECKPOINT_ROOT,
        "xtpijc3k",   # your run ID
        "epochs",
        "last.ckpt",
    )
    model = load_model(cfg, ckpt_path, gene_mapping)

    # ------------------------------------------------
    # RUN PREDICT()
    # ------------------------------------------------
    cls_emb, mean_emb = run_prediction(model, datamodule)

    # ------------------------------------------------
    # SAVE EMBEDDINGS — ONLY VALIDATION FILES
    # ------------------------------------------------
    val_files = cfg.PATH.SPLIT["val"]
    if isinstance(val_files, str):
        val_files = [val_files]

    for file in val_files:
        full_path = os.path.join(cfg.PATH.ADATA_PATH, file)
        save_embeddings(full_path, cls_emb, mean_emb)


if __name__ == "__main__":
    main()
