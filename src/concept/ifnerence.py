import os
import sys
import torch
import pandas as pd
import anndata as ad
import lightning as L

from concept.model import BiEncoderContrastiveModel
from concept.data.datamodules import AnnDataModule
from lamin_dataloader.dataset import GeneIdTokenizer
from omegaconf import OmegaConf
from hydra import compose, initialize


# ---------------------------
# LOAD CONFIG (same as train)
# ---------------------------
print("Loading config from YAML...")
with initialize(version_base=None, config_path="./conf"):
    cfg = compose(config_name="config", overrides=sys.argv[1:])

print(OmegaConf.to_yaml(cfg))

# ---------------------------
# LOAD CHECKPOINT + GENE MAP
# ---------------------------
print("\nLoading checkpoint and tokenizer...")

ckpt_path = os.path.join(
    cfg.PATH.CHECKPOINT_ROOT,
    "xtpijc3k",
    "epochs",
    "last.ckpt"
)

gene_mapping = pd.read_pickle(cfg.PATH.gene_mapping_path).to_dict()
tokenizer = GeneIdTokenizer(gene_mapping)

model_args = {
    "config": cfg.model,
    "pad_token_id": gene_mapping["<pad>"],
    "cls_token_id": gene_mapping["<cls>"],
    "vocab_size": len(gene_mapping),
    "world_size": 1,
    "val_loader_names": [],
    "precomp_embs_key": None,
}

print(f"Loading model from {ckpt_path}")
model = BiEncoderContrastiveModel.load_from_checkpoint(
    ckpt_path, **model_args, strict=False
)
model.eval()
model.cuda()


# ---------------------------
# BUILD DATAMODULE FOR NEW DATA
# ---------------------------
print("\nBuilding inference AnnDataModule...")

dataset_path = cfg.PATH.ADATA_PATH     # <--- YOUR .h5ad file


datamodule = AnnDataModule(
    dataset_path=dataset_path,
    split=cfg.PATH.SPLIT,              # use your new dataset
    panels_path=cfg.PATH.PANELS_PATH,
    columns=cfg.datamodule.columns,
    precomp_embs_key=None,
    normalization=cfg.datamodule.normalization,
    gene_sampling_strategy="none",              # IMPORTANT
    dataset_kwargs={"train": {"panel_selection": "preselected"}},
    dataloader_kwargs={"batch_size": 1, "num_workers": 4},
    val_loader_names=[],
    tokenizer=tokenizer,
)


# ---------------------------
# RUN PREDICTION
# ---------------------------
print("\nRunning model.predict()...")
trainer = L.Trainer(accelerator="gpu", devices=1, logger=False)
preds = trainer.predict(model, datamodule=datamodule)

print("Collecting embeddings...")
embeddings = torch.cat([p["cls_cell_emb"] for p in preds], dim=0).cpu().numpy()


# ---------------------------
# SAVE INTO ANNDATA
# ---------------------------
print("\nLoading AnnData and adding embedding...")
adata = ad.read_h5ad(adata_file)
adata.obsm["concept_embedding"] = embeddings

out_path = adata_file.replace(".h5ad", "_with_concept_embedding.h5ad")
adata.write(out_path)

print(f"\n🎉 Saved updated AnnData with embeddings to:\n{out_path}\n")
