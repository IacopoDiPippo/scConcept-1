"""
UMAP pipeline: loads embeddings like evaluate.py, builds adata_ref and plots like Code 2.

Usage:
    python umap_pipeline.py --model-dir /path/to/models/small_weighted__hgld4ax1_last \
                            --output-dir /path/to/output \
                            --subsample 50000
"""

import argparse
import os
import logging

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from typing import List, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


# ============================================================
# CONFIGURATION — edit paths here
# ============================================================

DATASET_CONFIGS = {
    "zeng": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zeng.h5ad",
        "annotation_path": "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zeng/results/merfish/cell_type_annotation/adata_obs_annotated.csv",
        "cell_type_col": "cell_type_mmc_raw",
        "cell_id_suffix": "-Zeng",         # suffix to strip from obs_names to match annotation
    },
    "zhuang": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zhuang-ABCA-1.h5ad",
        "annotation_path": "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv",
        "cell_type_col": "cell_type_mmc_raw",
        "cell_id_suffix": "-Zhuang-ABCA-1",
    },
    "isd": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/hvg_isd_normed.h5ad",
        "annotation_path": None,            # cell_type already in adata
        "cell_type_col": "cell_type_mmc_raw_revised",  # or "cell_type"
        "cell_id_suffix": None,
    },
}

CELL_TYPE_PALETTE = {
    'ABCs': '#023FA5', 'Astrocytes': '#7D87B9', 'Astroependymal': '#BEC1D4',
    'BAMs': '#D6BCC0', 'Bergmann': '#BB7784', 'Choroid-Plexus': '#8E063B',
    'ECs': '#4A6FE3', 'Ependymal': '#8595E1', 'Immune-Other': '#B5BBE3',
    'Microglia': '#E6AFB9', 'Neurons-Dopa': '#E07B91', 'Neurons-Gaba': '#D33F6A',
    'Neurons-Glut': '#11C638', 'Neurons-Glyc-Gaba': '#8DD593',
    'Neurons-Granule-Immature': '#C6DEC7', 'Neurons-Other': '#EAD3C6',
    'OECs': '#F0B98D', 'OPCs': '#EF9708', 'Oligodendrocytes': '#0FCFC0',
    'Pericytes': '#9CDED6', 'SMCs': '#D5EAE7', 'Tanycytes': '#F3E1EB',
    'Undefined': '#F6C4E1', 'Unknown': '#F6C4E1', 'VLMCs': '#F79CD4'
}


# ============================================================
# STEP 1: Load embeddings — exactly like evaluate.py
# ============================================================

def add_concept_embeddings(adata, embedding_path: str, name: str = "Dataset"):
    """
    Load concept embeddings from .npy files and align them to adata by cell_id.
    Exactly mirrors evaluate.py logic.
    """
    log.info(f"Loading embeddings for {name} from {embedding_path}")

    cell_ids = np.load(f"{embedding_path}/cell_ids.npy", allow_pickle=True).astype(str)
    emb_mean = np.load(f"{embedding_path}/cell_embs_mean.npy")
    emb_cls  = np.load(f"{embedding_path}/cell_embs_cls.npy")

    df_mean = pd.DataFrame(emb_mean, index=cell_ids)
    df_cls  = pd.DataFrame(emb_cls,  index=cell_ids)

    # Find cells present in both adata and embeddings
    common_cells = adata.obs_names.intersection(cell_ids)
    log.info(f"   {name}: {len(common_cells):,} / {adata.n_obs:,} cells have embeddings")

    adata = adata[common_cells].copy()

    # Align by obs_names — same as evaluate.py
    adata.obsm["concept_mean_embedding"] = df_mean.loc[adata.obs_names].to_numpy()
    adata.obsm["concept_cls_embedding"]  = df_cls.loc[adata.obs_names].to_numpy()

    return adata


# ============================================================
# STEP 2: Load dataset + annotations — like evaluate.py
# ============================================================

def load_dataset(name: str, model_dir: str, subsample: Optional[int] = None) -> ad.AnnData:
    """
    Load a single dataset, attach embeddings, attach cell_type annotations.
    """
    cfg = DATASET_CONFIGS[name]
    embedding_path = os.path.join(model_dir, name)

    log.info(f"\n--- Loading {name} ---")
    adata = sc.read_h5ad(cfg["adata_path"])

    # Drop X to save memory (like Code 2 does)
    adata.X = None

    # Add embeddings — this also filters to cells with embeddings
    adata = add_concept_embeddings(adata, embedding_path, name=name)

    # Add cell_type annotations
    if cfg["annotation_path"] is not None:
        ann = pd.read_csv(cfg["annotation_path"])

        # Strip suffix from obs_names to match annotation cell_ids
        suffix = cfg.get("cell_id_suffix")
        if suffix:
            obs_ids_stripped = adata.obs_names.str.replace(suffix, "", regex=False)
        else:
            obs_ids_stripped = adata.obs_names

        adata.obs["cell_id_stripped"] = obs_ids_stripped.values

        # Merge annotation
        ann = ann[["cell_id", cfg["cell_type_col"]]].rename(
            columns={cfg["cell_type_col"]: "cell_type", "cell_id": "cell_id_stripped"}
        )
        adata.obs = adata.obs.merge(ann, on="cell_id_stripped", how="left")
        adata.obs.index = adata.obs_names  # restore index after merge
    else:
        # ISD: cell_type already in adata, just rename to unified column
        col = cfg["cell_type_col"]
        if col in adata.obs.columns:
            adata.obs["cell_type"] = adata.obs[col]
        elif "cell_type" not in adata.obs.columns:
            log.warning(f"No cell_type column found for {name}")
            adata.obs["cell_type"] = "Unknown"

    # Standardise NaNs in cell_type
    adata.obs["cell_type"] = adata.obs["cell_type"].fillna("Unknown").astype(str).astype("category")

    # Add dataset label
    adata.obs["dataset"] = name.upper()

    # Subsample if requested
    if subsample and adata.n_obs > subsample:
        log.info(f"   Subsampling {name} from {adata.n_obs:,} to {subsample:,}")
        sc.pp.subsample(adata, n_obs=subsample, random_state=0)

    log.info(f"   Final: {adata.n_obs:,} cells")
    return adata


# ============================================================
# STEP 3: Compute UMAP and plot — exactly like Code 2
# ============================================================

def compute_and_plot_umap(
    adata: ad.AnnData,
    embedding_key: str,
    color_keys: List[str],
    output_dir: str,
    palette_map: Optional[dict] = None,
):
    log.info(f"Computing neighbors + UMAP using {embedding_key}...")
    sc.pp.neighbors(adata, n_neighbors=30, use_rep=embedding_key)
    sc.tl.umap(adata, min_dist=0.3, random_state=0)

    os.makedirs(output_dir, exist_ok=True)

    for color in color_keys:
        if color not in adata.obs.columns:
            log.warning(f"Skipping color '{color}': not in adata.obs")
            continue

        use_palette = None
        if palette_map:
            cats = adata.obs[color].unique()
            if any(c in palette_map for c in cats):
                use_palette = palette_map

        out_png = os.path.join(output_dir, f"umap_{embedding_key}_{color}.png")
        log.info(f"Saving: {out_png}")

        fig = sc.pl.umap(
            adata,
            color=color,
            title=f"UMAP {embedding_key} - {color}",
            palette=use_palette,
            frameon=False,
            show=False,
            return_fig=True,
        )
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        log.info(f"   Saved: {out_png}")


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="UMAP pipeline for Zeng/Zhuang/ISD")
    parser.add_argument(
        "--model-dir", required=True,
        help="Path to model embeddings directory, e.g. .../models/small_weighted__hgld4ax1_last"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Where to save UMAP plots"
    )
    parser.add_argument(
        "--subsample", type=int, default=50000,
        help="Cells per dataset after subsampling (default: 50000, set 0 to disable)"
    )
    parser.add_argument(
        "--embedding-key", default="concept_cls_embedding",
        choices=["concept_cls_embedding", "concept_mean_embedding"],
        help="Which embedding to use for UMAP"
    )
    parser.add_argument(
        "--datasets", nargs="+", default=["zeng", "zhuang", "isd"],
        choices=["zeng", "zhuang", "isd"],
        help="Datasets to include"
    )
    args = parser.parse_args()

    subsample = args.subsample if args.subsample > 0 else None

    # Step 1+2: Load all datasets
    adatas = {}
    for name in args.datasets:
        adatas[name] = load_dataset(name, args.model_dir, subsample=subsample)

    # Step 3: Concatenate — outer join, obs_names made unique (like evaluate.py)
    log.info("\nConcatenating datasets...")
    adata_ref = ad.concat(
        adatas,
        join="outer",
        label="dataset",
        index_unique="-",
    )
    adata_ref.obs_names_make_unique()

    log.info(f"Combined: {adata_ref.n_obs:,} cells")
    log.info(f"Dataset counts:\n{adata_ref.obs['dataset'].value_counts()}")

    # Step 4: UMAP + plots
    compute_and_plot_umap(
        adata=adata_ref,
        embedding_key=args.embedding_key,
        color_keys=["dataset", "cell_type"],
        output_dir=args.output_dir,
        palette_map=CELL_TYPE_PALETTE,
    )

    log.info("\nDone!")


if __name__ == "__main__":
    main()