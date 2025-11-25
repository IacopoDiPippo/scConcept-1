import os
import sys
import torch
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from concept.model import BiEncoderContrastiveModel
from lamin_dataloader.dataset import GeneIdTokenizer
from omegaconf import OmegaConf
from hydra import compose, initialize

# ---------------------------
# 1. Load config
# ---------------------------
print("Loading config from YAML...")
with initialize(version_base=None, config_path="./conf"):
    cfg = compose(config_name="config", overrides=sys.argv[1:])

print(OmegaConf.to_yaml(cfg))

# ---------------------------
# 2. Paths: checkpoint, gene mapping, adata file
# ---------------------------
print("\nLoading checkpoint and tokenizer...")

ckpt_path = os.path.join(
    cfg.PATH.CHECKPOINT_ROOT,
    "xtpijc3k",         # run id you gave earlier
    "epochs",
    "last.ckpt",
)

gene_mapping = pd.read_pickle(cfg.PATH.gene_mapping_path).to_dict()
tokenizer = GeneIdTokenizer(gene_mapping)

# IMPORTANT: this must be a SINGLE .h5ad file, not a directory
adata_file = cfg.PATH.ADATA_PATH
adata_file = os.path.join(adata_file, "Zeng.h5ad")  # in case of relative path
if os.path.isdir(adata_file):
    raise ValueError(
        f"cfg.PATH.ADATA_PATH points to a DIRECTORY:\n  {adata_file}\n"
        "You MUST set it to a specific .h5ad file, e.g.:\n"
        "  PATH.ADATA_PATH=/p/project1/.../Zeng.h5ad"
    )
if not adata_file.endswith(".h5ad"):
    raise ValueError(f"ADATA_PATH must be a .h5ad file, got: {adata_file}")
if not os.path.exists(adata_file):
    raise FileNotFoundError(f"File not found: {adata_file}")

print(f"Using AnnData file for inference:\n  {adata_file}")

# ---------------------------
# 3. Load model
# ---------------------------
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

dim_model = model.dim_model
print(f"Model embedding dim (dim_model): {dim_model}")

# ---------------------------
# 4. Load AnnData
# ---------------------------
print("\nLoading AnnData...")
adata = ad.read_h5ad(adata_file)
X = adata.X
var_names = np.array(adata.var_names)

print(f"AnnData shape: {adata.n_obs} cells × {adata.n_vars} genes")

# Map genes to token ids (vocabulary)
PAD_ID = gene_mapping["<pad>"]

gene_token_ids = np.full(adata.n_vars, PAD_ID, dtype=np.int64)
valid_gene_mask = np.zeros(adata.n_vars, dtype=bool)

for j, g in enumerate(var_names):
    if g in gene_mapping:
        gene_token_ids[j] = gene_mapping[g]
        valid_gene_mask[j] = True

print(f"Genes in vocab: {valid_gene_mask.sum()} / {adata.n_vars}")

# ---------------------------
# 5. Manual inference loop over cells
#    (no datamodule, no Trainer.predict)
# ---------------------------
max_tokens = cfg.datamodule.dataset.train.max_tokens  # 1000 in your config
print(f"Using up to {max_tokens} non-zero genes per cell.")

emb_list = []

for i in range(adata.n_obs):
    # Get expression vector
    if sparse.isspmatrix(X):
        x = X[i].toarray().ravel()
    else:
        x = np.asarray(X[i]).ravel()

    # Same normalization as in config (log1p)
    if cfg.datamodule.normalization == "log1p":
        x = np.log1p(x)

    # Select non-zero & valid genes
    mask = (x != 0) & valid_gene_mask
    idxs = np.where(mask)[0]

    if idxs.size == 0:
        # No usable genes for this cell -> use a dummy tiny embedding
        # (you may want to handle this differently)
        tokens = torch.tensor([[PAD_ID]], device=model.device, dtype=torch.long)
        values = torch.zeros_like(tokens, dtype=torch.float32)
    else:
        # Limit to max_tokens
        if idxs.size > max_tokens:
            idxs = np.random.choice(idxs, max_tokens, replace=False)

        tokens_1d = gene_token_ids[idxs]
        values_1d = x[idxs]

        tokens = torch.tensor(tokens_1d, device=model.device, dtype=torch.long).unsqueeze(0)
        values = torch.tensor(values_1d, device=model.device, dtype=torch.float32).unsqueeze(0)

    # No padding mask (mask_padding = False in your config)
    padding_mask = torch.zeros_like(tokens, dtype=torch.bool, device=model.device)

    with torch.no_grad():
        # forward: returns (pred, embs_padded, cell_embs)
        _, _, cell_emb = model(tokens, values, src_key_padding_mask=padding_mask)

    emb_list.append(cell_emb.cpu())

embeddings = torch.cat(emb_list, dim=0).numpy()
print("Final embedding matrix shape:", embeddings.shape)  # should be (n_cells, 128)

# ---------------------------
# 6. Save into AnnData in SCRATCH_ROOT
# ---------------------------
print("\nSaving embeddings into AnnData...")

adata.obsm["concept_embedding"] = embeddings

basename = os.path.basename(adata_file).replace(
    ".h5ad", "_with_concept_embedding.h5ad"
)
out_path = os.path.join(cfg.PATH.SCRATCH_ROOT, basename)

os.makedirs(cfg.PATH.SCRATCH_ROOT, exist_ok=True)
adata.write(out_path)

print(f"\n🎉 Saved updated AnnData with embeddings to:\n{out_path}\n")
