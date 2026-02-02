#!/usr/bin/env python3
"""
Pipeline for generating embeddings, computing UMAP, and calculating metrics (cLISI, iLISI)
for single-cell datasets (ISD, Zeng, Zhuang, Atlas, Atlas2, Atlas_train, Atlas2_train).

Usage examples (run from src/ directory):
    # Generate embeddings for all datasets
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets isd zeng zhuang

    # Generate embeddings + UMAP for specific combination
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets zeng zhuang --umap

    # Full pipeline with metrics
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets isd zeng zhuang \\
        --umap --clisi cell_type --ilisi dataset

    # Atlas with all panels + normal atlas
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets atlas \\
        --atlas-panel all --umap --ilisi dataset

    # Compare atlas and atlas2
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets atlas atlas2 \\
        --umap --ilisi dataset

    # Atlas_train with panels
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets atlas_train \\
        --atlas-train-panel zhuang zeng --umap --ilisi dataset

    # Atlas_train without panels (normal, all genes)
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets atlas_train \\
        --umap --ilisi dataset

    # Compare train vs val
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets atlas atlas_train \\
        --umap --ilisi dataset
"""

import argparse
from html import parser
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

# Set temp directory before importing scib
TMPDIR = "/p/scratch/cjinm16/dipippo1/tmp_lisi"
os.makedirs(TMPDIR, exist_ok=True)
os.environ["TMPDIR"] = TMPDIR
os.environ["TEMP"] = TMPDIR
os.environ["TMP"] = TMPDIR

import tempfile
tempfile.tempdir = TMPDIR


# ============================================================
# CONFIGURATION - Modify these paths as needed
# ============================================================

CONFIG = {
    # Dataset paths
    "isd": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/hvg_isd_normed.h5ad",
        "annotation_path": None,  # ISD has cell_type in adata already
        "has_cell_type": True,
        "subsample": 100_000,
    },
    "zeng": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zeng.h5ad",
        "annotation_path": "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zeng/results/merfish/cell_type_annotation/adata_obs_annotated.csv",
        "has_cell_type": True,
        "subsample": 100_000,
        "cell_id_suffix": "-Zeng",
    },
    "zhuang": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zhuang-ABCA-1.h5ad",
        "annotation_path": "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv",
        "has_cell_type": True,
        "subsample": 100_000,
        "cell_id_suffix": "-Zhuang-ABCA-1",
    },
    "atlas": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/train_val_data/abc_atlas_val.h5ad",
        "annotation_path": None,
        "has_cell_type": False,
        "subsample": 100_000,
    },
    "atlas2": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/train_val_data/abc_atlas_2_val.h5ad",
        "annotation_path": None,
        "has_cell_type": False,
        "subsample": 100_000,
    },
    "atlas_train": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/train_val_data/abc_atlas_train.h5ad",
        "annotation_path": None,
        "has_cell_type": False,
        "subsample": 100_000,
    },
    "atlas2_train": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/train_val_data/abc_atlas_2_train.h5ad",
        "annotation_path": None,
        "has_cell_type": False,
        "subsample": 100_000,
    },
}

# Panel paths for filtering atlas
PANELS = {
    "zhuang": "/p/scratch/cjinm16/dipippo1/scConcept/panels/Zhuang_ABCA1.csv",
    "zeng": "/p/scratch/cjinm16/dipippo1/scConcept/panels/Zeng_Panel.csv",
    "isd": "/p/scratch/cjinm16/dipippo1/scConcept/panels/ISD.csv",
}

# Path for filtered atlas files (separate directory)
FILTERED_ATLAS_PATH = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/atlas_filtered"

# Base path for embeddings (will be organized by model name)
EMBEDDINGS_BASE_PATH = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/models"

# Output figures path - UMAPs always saved here
FIGURES_PATH = "/p/scratch/cjinm16/dipippo1/scConcept/umaps"

# Cell type palette for visualization
CELL_TYPE_PALETTE = {
    "ABCs": "#023fa5",
    "Astrocytes": "#7d87b9",
    "Astroependymal": "#bec1d4",
    "BAMs": "#d6bcc0",
    "Bergmann": "#bb7784",
    "Choroid-Plexus": "#8e063b",
    "ECs": "#4a6fe3",
    "Ependymal": "#8595e1",
    "Immune-Other": "#b5bbe3",
    "Microglia": "#e6afb9",
    "Neurons-Dopa": "#e07b91",
    "Neurons-Gaba": "#d33f6a",
    "Neurons-Glut": "#11c638",
    "Neurons-Glyc-Gaba": "#8dd593",
    "Neurons-Granule-Immature": "#c6dec7",
    "Neurons-Other": "#ead3c6",
    "OECs": "#f0b98d",
    "OPCs": "#ef9708",
    "Oligodendrocytes": "#0fcfc0",
    "Pericytes": "#9cded6",
    "SMCs": "#d5eae7",
    "Tanycytes": "#f3e1eb",
    "Undefined": "#f6c4e1",
    "VLMCs": "#f79cd4",
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_model_name(checkpoint_path: str) -> str:
    """Extract model name from checkpoint path for organizing embeddings."""
    path = Path(checkpoint_path)
    step_name = path.stem.replace("=", "")  # step=310000 -> step310000
    model_dir = path.parent.parent.name  # zp2ksa3s
    return f"{model_dir}_{step_name}"


def get_embedding_path(model_name: str, dataset_name: str, panel: str = None) -> str:
    """Get the embedding output path for a given model and dataset."""
    if panel:
        return os.path.join(EMBEDDINGS_BASE_PATH, model_name, f"{dataset_name}_{panel}")
    return os.path.join(EMBEDDINGS_BASE_PATH, model_name, dataset_name)


def get_filtered_atlas_path(atlas_name: str, panel: str) -> str:
    """Get path for panel-filtered atlas file."""
    base_name = get_atlas_base_name(atlas_name)
    return os.path.join(FILTERED_ATLAS_PATH, f"{base_name}_{panel}.h5ad")


def filter_adata_by_panel(adata_path: str, panel_path: str, output_path: str) -> str:
    """Filter AnnData object to only include genes from a panel."""
    import anndata as ad
    
    if os.path.exists(output_path):
        print(f"✅ Filtered atlas already exists: {output_path}")
        return output_path
    
    print(f"🔧 Filtering atlas by panel...")
    print(f"   Input: {adata_path}")
    print(f"   Panel: {panel_path}")
    
    adata = ad.read_h5ad(adata_path)
    panel = pd.read_csv(panel_path)
    
    if 'Ensembl_ID' not in panel.columns:
        raise ValueError(f"Panel CSV must have 'Ensembl_ID' column. Found: {panel.columns.tolist()}")
    
    panel_genes = panel['Ensembl_ID'].tolist()
    
    # Filter to keep only genes in the panel
    adata_filtered = adata[:, adata.var_names.isin(panel_genes)]
    genes_found = adata.var_names.isin(panel_genes).sum()
    
    print(f"   Original: {adata.shape[1]} genes -> Filtered: {adata_filtered.shape[1]} genes")
    print(f"   Found {genes_found}/{len(panel_genes)} panel genes")
    
    adata_filtered.write_h5ad(output_path)
    print(f"   ✅ Saved: {output_path}")
    
    return output_path


def embeddings_exist(embedding_path: str) -> bool:
    """Check if embeddings already exist."""
    required_files = ["cell_ids.npy", "cell_embs_mean.npy", "cell_embs_cls.npy"]
    return all(os.path.exists(os.path.join(embedding_path, f)) for f in required_files)


def generate_embeddings(checkpoint: str, adata_path: str, output_path: str, batch_size: int = 64, subsample: int = None):
    """Generate embeddings using concept.get_embs module."""
    print(f"\n{'='*60}")
    print(f"🚀 Generating embeddings")
    print(f"   Checkpoint: {checkpoint}")
    print(f"   Input: {adata_path}")
    print(f"   Output: {output_path}")
    if subsample:
        print(f"   Subsample: {subsample:,} cells")
    print(f"{'='*60}\n")
    
    os.makedirs(output_path, exist_ok=True)
    
    cmd = [
        "python", "-m", "concept.get_embs",
        "--checkpoint", checkpoint,
        "--adata_path", adata_path,
        "--output_emb_path", output_path,
        "--batch_size", str(batch_size),
    ]
    
    # AGGIUNGI SUBSAMPLE SE SPECIFICATO
    if subsample:
        cmd.extend(["--subsample", str(subsample)])
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        raise RuntimeError(f"Embedding generation failed with return code {result.returncode}")
    
    print(f"✅ Embeddings saved to: {output_path}")


def sanity_check(adata, name="adata"):
    """Print minimal info about an AnnData object."""
    print(f"📋 {name}: {adata.n_obs:,} cells × {adata.n_vars:,} genes | obsm: {list(adata.obsm.keys())}")


def add_concept_embeddings(adata, embedding_path: str, name: str = "Dataset"):
    """Load concept embeddings and add them to adata.obsm."""
    cell_ids = np.load(f"{embedding_path}/cell_ids.npy", allow_pickle=True).astype(str)
    emb_mean = np.load(f"{embedding_path}/cell_embs_mean.npy")
    emb_cls = np.load(f"{embedding_path}/cell_embs_cls.npy")
    
    df_mean = pd.DataFrame(emb_mean, index=cell_ids)
    df_cls = pd.DataFrame(emb_cls, index=cell_ids)
    
    # Align with adata.obs_names
    adata.obsm["concept_mean_embedding"] = df_mean.loc[adata.obs_names].to_numpy()
    adata.obsm["concept_cls_embedding"] = df_cls.loc[adata.obs_names].to_numpy()
    
    return adata


def load_dataset(dataset_name: str, embedding_path: str, subsample: Optional[int] = None):
    """Load a dataset with its embeddings and annotations."""
    import anndata as ad
    import scanpy as sc
    
    config = CONFIG[dataset_name]
    print(f"📂 Loading {dataset_name}...")
    
    adata = sc.read(config["adata_path"])
    
    # Add embeddings
    adata = add_concept_embeddings(adata, embedding_path, name=dataset_name)
    
    # Handle annotations if available
    if config["annotation_path"] is not None:
        ann = pd.read_csv(config["annotation_path"])
        
        # Create cell_id with proper suffix for adata
        suffix = config.get("cell_id_suffix", f"-{dataset_name}")
        adata.obs["cell_id"] = adata.obs.index.astype(str) + suffix
        
        # For Zeng, also add suffix to annotation CSV (it doesn't have it)
        # For Zhuang, the annotation CSV already has the full cell_id format
        if dataset_name == "zeng":
            ann["cell_id"] = ann["cell_id"].astype(str) + suffix
        
        # Merge annotations
        adata.obs = adata.obs.merge(
            ann[["cell_id", "cell_type_mmc_raw"]],
            on="cell_id",
            how="left"
        )
        
        # Filter to annotated cells only
        mask_valid = adata.obs["cell_type_mmc_raw"].notna()
        print(f"   Annotated: {mask_valid.sum()}/{adata.n_obs} cells")
        adata = adata[mask_valid].copy()
    
    # For ISD, check if cell_type column exists
    if dataset_name == "isd":
        if "cell_type_mmc_raw_revised" in adata.obs.columns:
            adata.obs["cell_type_mmc_raw"] = adata.obs["cell_type_mmc_raw_revised"]
        elif "cell_type" in adata.obs.columns:
            adata.obs["cell_type_mmc_raw"] = adata.obs["cell_type"]
    
    # Subsample if needed
    if subsample and adata.n_obs > subsample:
        print(f"   Subsampling to {subsample} cells")
        sc.pp.subsample(adata, n_obs=subsample, random_state=0)
    
    # Add dataset identifier
    adata.obs["dataset"] = dataset_name.capitalize()
    
    print(f"   ✅ {adata.n_obs:,} cells loaded")
    return adata


def load_dataset_custom(dataset_name: str, adata_path: str, embedding_path: str, 
                        subsample: Optional[int] = None, has_cell_type: bool = False):
    """Load a dataset from custom path (for filtered atlas)."""
    import anndata as ad
    import scanpy as sc
    
    print(f"📂 Loading {dataset_name} (custom)...")
    
    adata = sc.read(adata_path)
    
    # Add embeddings
    adata = add_concept_embeddings(adata, embedding_path, name=dataset_name)
    
    # Subsample if needed
    if subsample and adata.n_obs > subsample:
        print(f"   Subsampling to {subsample} cells")
        sc.pp.subsample(adata, n_obs=subsample, random_state=0)
    
    # Add dataset identifier
    adata.obs["dataset"] = dataset_name.capitalize()
    
    print(f"   ✅ {adata.n_obs:,} cells loaded")
    return adata


def compute_neighbors_umap(adata, use_rep: str = "concept_cls_embedding"):
    """Compute neighbors and UMAP."""
    import scanpy as sc
    
    print(f"🔄 Computing neighbors using {use_rep}...")
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=30)
    
    print("🔄 Computing UMAP...")
    sc.tl.umap(adata, min_dist=0.3, random_state=0)


def plot_umap(adata, color: str, palette=None, title: str = None, out_png: str = None):
    """Plot UMAP and save to file."""
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    fig = sc.pl.umap(
        adata,
        color=color,
        title=title,
        frameon=False,
        palette=palette,
        show=False,
        return_fig=True,
    )
    
    if out_png:
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        print(f"💾 Saved: {out_png}")
    
    plt.close(fig)


def compute_clisi(adata, label_key: str) -> float:
    """Compute cLISI score."""
    import scib
    
    print(f"📊 Computing cLISI for '{label_key}'...")
    score = scib.metrics.clisi_graph(adata, label_key=label_key, type_="knn")
    print(f"   cLISI ({label_key}): {score:.4f}")
    return score


def compute_ilisi(adata, batch_key: str) -> float:
    """Compute iLISI score."""
    import scib
    
    print(f"📊 Computing iLISI for '{batch_key}'...")
    score = scib.metrics.ilisi_graph(adata, batch_key=batch_key, type_="knn")
    print(f"   iLISI ({batch_key}): {score:.4f}")
    return score


def is_atlas_dataset(ds: str) -> bool:
    """Check if dataset is an atlas type."""
    return ds in ["atlas", "atlas2", "atlas_train", "atlas2_train"]


def get_atlas_base_name(ds: str) -> str:
    """Get the base atlas name for filtered file naming."""
    mapping = {
        "atlas": "abc_atlas_val",
        "atlas2": "abc_atlas_2_val",
        "atlas_train": "abc_atlas_train",
        "atlas2_train": "abc_atlas_2_train",
    }
    return mapping.get(ds, ds)


# ============================================================
# MAIN PIPELINE
# ============================================================

def run_pipeline(
    checkpoint: str,
    datasets: List[str],
    compute_umap: bool = False,
    clisi_keys: List[str] = None,
    ilisi_keys: List[str] = None,
    batch_size: int = 64,
    force_regenerate: bool = False,
    atlas_panels: List[str] = None,
    atlas2_panels: List[str] = None,
    atlas_train_panels: List[str] = None,
    atlas2_train_panels: List[str] = None,
):
    """Run the full pipeline."""
    import anndata as ad
    import scanpy as sc
    
    model_name = get_model_name(checkpoint)
    
    # Map dataset to its panels argument
    panels_map = {
        "atlas": atlas_panels,
        "atlas2": atlas2_panels,
        "atlas_train": atlas_train_panels,
        "atlas2_train": atlas2_train_panels,
    }
    
    # Handle 'all' option for each atlas type and track which need normal version
    include_normal = {}
    for ds_name in ["atlas", "atlas2", "atlas_train", "atlas2_train"]:
        panels = panels_map.get(ds_name)
        include_normal[ds_name] = False
        if panels and "all" in panels:
            panels_map[ds_name] = ["zhuang", "zeng", "isd"]
            include_normal[ds_name] = True
    
    # Update local variables
    atlas_panels = panels_map["atlas"]
    atlas2_panels = panels_map["atlas2"]
    atlas_train_panels = panels_map["atlas_train"]
    atlas2_train_panels = panels_map["atlas2_train"]
    
    # Build panel info string
    panel_info = ""
    for ds_name, panels in panels_map.items():
        if panels and ds_name in datasets:
            panel_info += f" ({ds_name} panels: {', '.join(panels)})"
            if include_normal[ds_name]:
                panel_info += " + normal"
    
    print(f"\n🚀 Pipeline: {model_name} | Datasets: {', '.join(datasets)}{panel_info}\n")
    
    # Validate datasets
    for ds in datasets:
        if ds not in CONFIG:
            raise ValueError(f"Unknown dataset: {ds}. Available: {list(CONFIG.keys())}")
    
    # Validate panels - warn if panel specified but dataset not in list
    for ds_name, panels in panels_map.items():
        if panels and ds_name not in datasets:
            print(f"⚠️ Warning: --{ds_name.replace('_', '-')}-panel specified but '{ds_name}' not in datasets. Ignoring.")
            panels_map[ds_name] = None
            include_normal[ds_name] = False
    
    # Re-update local variables after validation
    atlas_panels = panels_map["atlas"]
    atlas2_panels = panels_map["atlas2"]
    atlas_train_panels = panels_map["atlas_train"]
    atlas2_train_panels = panels_map["atlas2_train"]
    
    # Validate panel names
    for ds_name, panels in panels_map.items():
        if panels:
            for panel in panels:
                if panel not in PANELS:
                    raise ValueError(f"Unknown panel: {panel}. Available: {list(PANELS.keys())}")
    
    # Helper to get panels for a dataset
    def get_panels_for_ds(ds):
        return panels_map.get(ds)
    
    def get_include_normal(ds):
        return include_normal.get(ds, False)
    
    # Step 1: Generate embeddings if needed
    print("=" * 60)
    print("STEP 1: EMBEDDINGS")
    print("=" * 60)
    
    for ds in datasets:
        ds_panels = get_panels_for_ds(ds)
        
        if is_atlas_dataset(ds) and ds_panels:
            # Generate embeddings for each panel
            for panel in ds_panels:
                emb_path = get_embedding_path(model_name, ds, panel=panel)
                
                if embeddings_exist(emb_path) and not force_regenerate:
                    print(f"✅ {ds}_{panel}: embeddings exist at {emb_path}")
                else:
                    # First, create filtered atlas if needed
                    filtered_atlas_path = get_filtered_atlas_path(ds, panel)
                    os.makedirs(os.path.dirname(filtered_atlas_path), exist_ok=True)
                    filter_adata_by_panel(
                        adata_path=CONFIG[ds]["adata_path"],
                        panel_path=PANELS[panel],
                        output_path=filtered_atlas_path,
                    )
                    # Then generate embeddings
                    generate_embeddings(
                        checkpoint=checkpoint,
                        adata_path=filtered_atlas_path,
                        output_path=emb_path,
                        batch_size=batch_size,
                        subsample=CONFIG[ds].get("subsample"),
                    )
            
            # Also generate normal atlas if requested
            if get_include_normal(ds):
                emb_path = get_embedding_path(model_name, ds)
                if embeddings_exist(emb_path) and not force_regenerate:
                    print(f"✅ {ds}: embeddings exist at {emb_path}")
                else:
                    generate_embeddings(
                        checkpoint=checkpoint,
                        adata_path=CONFIG[ds]["adata_path"],
                        output_path=emb_path,
                        batch_size=batch_size,
                        subsample=CONFIG[ds].get("subsample"),
                    )
        
        else:
            # Normal dataset or atlas without panels
            emb_path = get_embedding_path(model_name, ds)
            
            if embeddings_exist(emb_path) and not force_regenerate:
                print(f"✅ {ds}: embeddings exist at {emb_path}")
            else:
                generate_embeddings(
                    checkpoint=checkpoint,
                    adata_path=CONFIG[ds]["adata_path"],
                    output_path=emb_path,
                    batch_size=batch_size,
                    subsample=CONFIG[ds].get("subsample"),
                )
    
    # Step 2: Load datasets
    print("\n" + "=" * 60)
    print("STEP 2: LOAD DATASETS")
    print("=" * 60)
    
    adatas = {}
    for ds in datasets:
        ds_panels = get_panels_for_ds(ds)
        
        if is_atlas_dataset(ds) and ds_panels:
            # Load each panel-filtered atlas separately
            for panel in ds_panels:
                emb_path = get_embedding_path(model_name, ds, panel=panel)
                filtered_atlas_path = get_filtered_atlas_path(ds, panel)
                key = f"{ds}_{panel}"
                adatas[key] = load_dataset_custom(
                    dataset_name=key,
                    adata_path=filtered_atlas_path,
                    embedding_path=emb_path,
                    subsample=CONFIG[ds].get("subsample"),
                    has_cell_type=False,
                )
            
            # Also load normal atlas if requested
            if get_include_normal(ds):
                emb_path = get_embedding_path(model_name, ds)
                adatas[ds] = load_dataset_custom(
                    dataset_name=ds,
                    adata_path=CONFIG[ds]["adata_path"],
                    embedding_path=emb_path,
                    subsample=CONFIG[ds].get("subsample"),
                    has_cell_type=False,
                )
        
        elif not is_atlas_dataset(ds):
            # Normal dataset (zeng, zhuang, isd)
            emb_path = get_embedding_path(model_name, ds)
            adatas[ds] = load_dataset(
                dataset_name=ds,
                embedding_path=emb_path,
                subsample=CONFIG[ds].get("subsample"),
            )
        
        else:
            # Atlas without panels
            emb_path = get_embedding_path(model_name, ds)
            adatas[ds] = load_dataset_custom(
                dataset_name=ds,
                adata_path=CONFIG[ds]["adata_path"],
                embedding_path=emb_path,
                subsample=CONFIG[ds].get("subsample"),
                has_cell_type=False,
            )
    
    # Step 3: Combine datasets
    print("\n" + "=" * 60)
    print("STEP 3: COMBINE")
    print("=" * 60)
    
    if len(adatas) == 1:
        combined = list(adatas.values())[0]
    else:
        combined = ad.concat(
            adatas,
            join="outer",
            label="dataset",
            index_unique=None,
        )
        combined.obs_names_make_unique()
    
    print(f"Combined: {combined.n_obs:,} cells | {dict(combined.obs['dataset'].value_counts())}")
    
    # Step 4: Compute neighbors (needed for UMAP and LISI metrics)
    if compute_umap or clisi_keys or ilisi_keys:
        print("\n" + "=" * 60)
        print("STEP 4: NEIGHBORS + UMAP")
        print("=" * 60)
        
        sc.pp.neighbors(combined, use_rep="concept_cls_embedding", n_neighbors=30)
    
    # Step 5: UMAP - always save to scratch
    if compute_umap:
        sc.tl.umap(combined, min_dist=0.3, random_state=0)
        
        # Build filename
        umap_dir = os.path.join(FIGURES_PATH, model_name)
        
        # Build dataset string for filename
        datasets_parts = []
        for ds in datasets:
            ds_panels = get_panels_for_ds(ds)
            if is_atlas_dataset(ds) and ds_panels:
                for panel in ds_panels:
                    datasets_parts.append(f"{ds}-{panel}")
                if get_include_normal(ds):
                    datasets_parts.append(ds)
            else:
                datasets_parts.append(ds)
        datasets_str = "_".join(datasets_parts)
        
        # Plot by dataset
        plot_umap(
            combined,
            color="dataset",
            title=f"UMAP by Dataset ({model_name})",
            out_png=os.path.join(umap_dir, f"{datasets_str}_dataset.png"),
        )
        
        # Plot by cell_type only if we have datasets with cell_type (not atlas-only)
        non_atlas_datasets = [ds for ds in datasets if not is_atlas_dataset(ds)]
        has_cell_type = any(CONFIG[ds]["has_cell_type"] for ds in non_atlas_datasets) if non_atlas_datasets else False
        if has_cell_type and "cell_type_mmc_raw" in combined.obs.columns:
            plot_umap(
                combined,
                color="cell_type_mmc_raw",
                palette=CELL_TYPE_PALETTE,
                title=f"UMAP by Cell Type ({model_name})",
                out_png=os.path.join(umap_dir, f"{datasets_str}_celltype.png"),
            )
    
    # Step 6: Compute metrics
    results = {"model": model_name, "datasets": "_".join(datasets)}
    for ds_name, panels in panels_map.items():
        if panels:
            results[f"{ds_name}_panels"] = ",".join(panels)
    
    if clisi_keys:
        print("\n" + "=" * 60)
        print("STEP 5: cLISI")
        print("=" * 60)
        
        for key in clisi_keys:
            col = "dataset" if key == "dataset" else "cell_type_mmc_raw"
            if col not in combined.obs.columns:
                print(f"⚠️ Skipping cLISI for {key}: column not found")
                continue
            results[f"cLISI_{key}"] = compute_clisi(combined, col)
    
    if ilisi_keys:
        print("\n" + "=" * 60)
        print("STEP 6: iLISI")
        print("=" * 60)
        
        for key in ilisi_keys:
            col = "dataset" if key == "dataset" else "cell_type_mmc_raw"
            if col not in combined.obs.columns:
                print(f"⚠️ Skipping iLISI for {key}: column not found")
                continue
            results[f"iLISI_{key}"] = compute_ilisi(combined, col)
    
    # Summary
    print("\n" + "=" * 60)
    print("📋 RESULTS")
    print("=" * 60)
    for k, v in results.items():
        print(f"   {k}: {v:.4f}" if isinstance(v, float) else f"   {k}: {v}")
    
    return combined, results


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Embedding generation and analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    
    parser.add_argument(
        "--checkpoint", "-c",
        required=True,
        help="Path to model checkpoint"
    )
    
    parser.add_argument(
        "--datasets", "-d",
        nargs="+",
        choices=["isd", "zeng", "zhuang", "atlas", "atlas2", "atlas_train", "atlas2_train"],
        required=True,
        help="Datasets to process"
    )
    
    parser.add_argument(
        "--umap",
        action="store_true",
        help="Compute and plot UMAP"
    )
    
    parser.add_argument(
        "--clisi",
        nargs="*",
        choices=["dataset", "cell_type"],
        help="Compute cLISI for specified keys"
    )
    
    parser.add_argument(
        "--ilisi",
        nargs="*",
        choices=["dataset", "cell_type"],
        help="Compute iLISI for specified keys"
    )
    
    parser.add_argument(
        "--batch-size",
        type=int,
        default=64,
        help="Batch size for embedding generation"
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force regeneration of embeddings even if they exist"
    )
    
    parser.add_argument(
        "--atlas-panel",
        nargs="+",
        choices=["zhuang", "zeng", "isd", "all"],
        help="Filter atlas to specific gene panel(s). Use 'all' to include atlas_zhuang, atlas_zeng, atlas_isd AND atlas normal"
    )
    
    parser.add_argument(
        "--atlas2-panel",
        nargs="+",
        choices=["zhuang", "zeng", "isd", "all"],
        help="Filter atlas2 to specific gene panel(s). Use 'all' to include atlas2_zhuang, atlas2_zeng, atlas2_isd AND atlas2 normal"
    )
    
    parser.add_argument(
        "--atlas-train-panel",
        nargs="+",
        choices=["zhuang", "zeng", "isd", "all"],
        help="Filter atlas_train to specific gene panel(s). Use 'all' for all panels + normal"
    )
    
    parser.add_argument(
        "--atlas2-train-panel",
        nargs="+",
        choices=["zhuang", "zeng", "isd", "all"],
        help="Filter atlas2_train to specific gene panel(s). Use 'all' for all panels + normal"
    )

    # Aggiungi questo argomento nel parser
    parser.add_argument("--subsample", type=int, default=None, help="Subsample to N cells before embedding")
    
    args = parser.parse_args()
    
    run_pipeline(
        checkpoint=args.checkpoint,
        datasets=args.datasets,
        compute_umap=args.umap,
        clisi_keys=args.clisi,
        ilisi_keys=args.ilisi,
        batch_size=args.batch_size,
        force_regenerate=args.force,
        atlas_panels=args.atlas_panel,
        atlas2_panels=args.atlas2_panel,
        atlas_train_panels=args.atlas_train_panel,
        atlas2_train_panels=args.atlas2_train_panel,
    )


if __name__ == "__main__":
    main()