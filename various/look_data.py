#!/usr/bin/env python3
import scanpy as sc
import anndata as ad
from pathlib import Path

# === Paths ===
DATA_DIR = Path("/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw")
FILES = ["abc_atlas.h5ad", "Zeng.h5ad"]

def summarize_adata(adata: ad.AnnData, name: str):
    print(f"\n📘 === Summary for {name} ===")
    print(f"Shape: {adata.shape}")
    print(f"n_obs (cells): {adata.n_obs}")
    print(f"n_vars (genes): {adata.n_vars}")

    # Basic keys
    print(f"Layers: {list(adata.layers.keys()) if hasattr(adata, 'layers') else 'None'}")
    print(f"Obs columns: {list(adata.obs.columns)}")
    print(f"Var columns: {list(adata.var.columns)}")
    print(f"Obsm keys: {list(adata.obsm.keys())}")
    print(f"Uns keys: {list(adata.uns.keys())}")

    # Small previews
    print("\n.obs head:")
    print(adata.obs.head())

    print("\n.var head:")
    print(adata.var.head())

def main():
    for fname in FILES:
        fpath = DATA_DIR / fname
        if not fpath.exists():
            print(f"❌ File not found: {fpath}")
            continue

        print(f"\n📂 Loading {fpath} (read-only)...")
        # Using mode='r' ensures read-only access, avoids copying data to memory unnecessarily
        adata = ad.read_h5ad(fpath, backed='r')

        summarize_adata(adata, fname)

if __name__ == "__main__":
    main()
