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

    # Atlas with ALL panels from directory (no normal)
    python -m concept.evaluate --checkpoint /path/to/model.ckpt --datasets atlas \\
        --atlas-panel all_all -s 50000 --umap

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

    # Force regenerate embeddings
    python -m concept.evaluate -c /path/to/checkpoint.ckpt -d zeng zhuang isd \\
        -s 100000 --umap --clisi cell_type --ilisi dataset --force
"""

import argparse
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
    "isd": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/hvg_isd_normed.h5ad",
        "annotation_path": None,
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

# Panel paths - load dynamically from directory
PANELS_DIR = "/p/scratch/cjinm16/dipippo1/scConcept/panels"

def load_panels_from_dir(panels_dir: str) -> dict:
    """Load all panel CSV files from directory."""
    panels = {}
    if os.path.exists(panels_dir):
        for f in os.listdir(panels_dir):
            if f.endswith('.csv'):
                panel_name = f.replace('.csv', '').lower().split('_')[0]
                panels[panel_name] = os.path.join(panels_dir, f)
    return panels

PANELS = load_panels_from_dir(PANELS_DIR)

# Path for filtered atlas files
FILTERED_ATLAS_PATH = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/atlas_filtered"

# Base path for embeddings
EMBEDDINGS_BASE_PATH = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/models"

# Output figures path
FIGURES_PATH = "/p/scratch/cjinm16/dipippo1/scConcept/umaps"

# Cell type palette
CELL_TYPE_PALETTE = {
    "ABCs": "#023fa5", "Astrocytes": "#7d87b9", "Astroependymal": "#bec1d4",
    "BAMs": "#d6bcc0", "Bergmann": "#bb7784", "Choroid-Plexus": "#8e063b",
    "ECs": "#4a6fe3", "Ependymal": "#8595e1", "Immune-Other": "#b5bbe3",
    "Microglia": "#e6afb9", "Neurons-Dopa": "#e07b91", "Neurons-Gaba": "#d33f6a",
    "Neurons-Glut": "#11c638", "Neurons-Glyc-Gaba": "#8dd593",
    "Neurons-Granule-Immature": "#c6dec7", "Neurons-Other": "#ead3c6",
    "OECs": "#f0b98d", "OPCs": "#ef9708", "Oligodendrocytes": "#0fcfc0",
    "Pericytes": "#9cded6", "SMCs": "#d5eae7", "Tanycytes": "#f3e1eb",
    "Undefined": "#f6c4e1", "VLMCs": "#f79cd4",
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_checkpoint_dir(checkpoint_path: str) -> Path:
    """Get checkpoint directory containing run_info.yaml and overrides.txt."""
    path = Path(checkpoint_path)
    if path.parent.name in ['steps', 'epochs']:
        return path.parent.parent
    return path.parent


def get_hydra_overrides(checkpoint_path: str) -> List[str]:
    """Read overrides.txt from checkpoint directory."""
    ckpt_dir = get_checkpoint_dir(checkpoint_path)
    overrides_file = ckpt_dir / "overrides.txt"
    
    if not overrides_file.exists():
        print(f"⚠️ WARNING: No overrides.txt found at {overrides_file}")
        return []
    
    with open(overrides_file, 'r') as f:
        overrides_str = f.read().strip()
    
    overrides = overrides_str.split()
    
    # Filter out overrides not needed for inference
    filtered = []
    skip_prefixes = (
        'wandb.',
        '+initialize_2.',
        'datamodule.dataloader.',  # batch_size etc not needed for inference
        'datamodule.probabilistic_panel_sampling',  # not needed for inference
    )
    for o in overrides:
        if any(o.startswith(prefix) for prefix in skip_prefixes):
            continue
        filtered.append(o)
    
    print(f"📋 Loaded overrides from {overrides_file}:")
    for o in filtered:
        print(f"   {o}")
    
    return filtered


def get_model_name(checkpoint_path: str) -> str:
    """Extract model name from checkpoint path."""
    path = Path(checkpoint_path)
    step_name = path.stem.replace("=", "")
    ckpt_dir = get_checkpoint_dir(checkpoint_path)
    return f"{ckpt_dir.name}_{step_name}"


def get_embedding_path(model_name: str, dataset_name: str, panel: str = None) -> str:
    """Get embedding output path."""
    if panel:
        return os.path.join(EMBEDDINGS_BASE_PATH, model_name, f"{dataset_name}_{panel}")
    return os.path.join(EMBEDDINGS_BASE_PATH, model_name, dataset_name)


def get_filtered_atlas_path(atlas_name: str, panel: str) -> str:
    """Get path for panel-filtered atlas file."""
    base_name = get_atlas_base_name(atlas_name)
    return os.path.join(FILTERED_ATLAS_PATH, f"{base_name}_{panel}.h5ad")


def filter_adata_by_panel(adata_path: str, panel_path: str, output_path: str) -> str:
    """Filter AnnData to include only genes from a panel."""
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


def generate_embeddings(
    checkpoint: str,
    adata_path: str,
    output_path: str,
    hydra_overrides: List[str],
    batch_size: int = 64,
    subsample: int = None
):
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
    
    if subsample:
        cmd.extend(["--subsample", str(subsample)])
    
    cmd.extend(hydra_overrides)
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        raise RuntimeError(f"Embedding generation failed with return code {result.returncode}")
    
    print(f"✅ Embeddings saved to: {output_path}")


def add_concept_embeddings(adata, embedding_path: str, name: str = "Dataset"):
    """Load concept embeddings and add them to adata.obsm."""
    cell_ids = np.load(f"{embedding_path}/cell_ids.npy", allow_pickle=True).astype(str)
    emb_mean = np.load(f"{embedding_path}/cell_embs_mean.npy")
    emb_cls = np.load(f"{embedding_path}/cell_embs_cls.npy")
    
    df_mean = pd.DataFrame(emb_mean, index=cell_ids)
    df_cls = pd.DataFrame(emb_cls, index=cell_ids)
    
    common_cells = adata.obs_names.intersection(cell_ids)
    print(f"   {name}: {len(common_cells):,}/{adata.n_obs:,} cells have embeddings")
    adata = adata[common_cells].copy()
    
    adata.obsm["concept_mean_embedding"] = df_mean.loc[adata.obs_names].to_numpy()
    adata.obsm["concept_cls_embedding"] = df_cls.loc[adata.obs_names].to_numpy()
    
    return adata


def load_dataset(dataset_name: str, embedding_path: str, subsample: Optional[int] = None):
    """Load a dataset with embeddings and annotations."""
    import scanpy as sc
    
    config = CONFIG[dataset_name]
    print(f"📂 Loading {dataset_name}...")
    
    adata = sc.read(config["adata_path"])
    adata = add_concept_embeddings(adata, embedding_path, name=dataset_name)
    
    if config["annotation_path"] is not None:
        ann = pd.read_csv(config["annotation_path"])
        suffix = config.get("cell_id_suffix", f"-{dataset_name}")
        adata.obs["cell_id"] = adata.obs.index.astype(str) + suffix
        
        if dataset_name == "zeng":
            ann["cell_id"] = ann["cell_id"].astype(str) + suffix
        
        adata.obs = adata.obs.merge(ann[["cell_id", "cell_type_mmc_raw"]], on="cell_id", how="left")
        mask_valid = adata.obs["cell_type_mmc_raw"].notna()
        print(f"   Annotated: {mask_valid.sum()}/{adata.n_obs} cells")
        adata = adata[mask_valid].copy()
    
    if dataset_name == "isd":
        if "cell_type_mmc_raw_revised" in adata.obs.columns:
            adata.obs["cell_type_mmc_raw"] = adata.obs["cell_type_mmc_raw_revised"]
        elif "cell_type" in adata.obs.columns:
            adata.obs["cell_type_mmc_raw"] = adata.obs["cell_type"]
    
    if subsample and adata.n_obs > subsample:
        print(f"   Subsampling to {subsample} cells")
        sc.pp.subsample(adata, n_obs=subsample, random_state=0)
    
    adata.obs["dataset"] = dataset_name.capitalize()
    print(f"   ✅ {adata.n_obs:,} cells loaded")
    return adata


def load_dataset_custom(dataset_name: str, adata_path: str, embedding_path: str,
                        subsample: Optional[int] = None, has_cell_type: bool = False):
    """Load a dataset from custom path (for filtered atlas)."""
    import scanpy as sc
    
    print(f"📂 Loading {dataset_name} (custom)...")
    adata = sc.read(adata_path)
    adata = add_concept_embeddings(adata, embedding_path, name=dataset_name)
    
    if subsample and adata.n_obs > subsample:
        print(f"   Subsampling to {subsample} cells")
        sc.pp.subsample(adata, n_obs=subsample, random_state=0)
    
    adata.obs["dataset"] = dataset_name.capitalize()
    print(f"   ✅ {adata.n_obs:,} cells loaded")
    return adata


def plot_umap(adata, color: str, palette=None, title: str = None, out_png: str = None):
    """Plot UMAP and save to file."""
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    fig = sc.pl.umap(adata, color=color, title=title, frameon=False,
                     palette=palette, show=False, return_fig=True)
    
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
    subsample: Optional[int] = None,
    mean: bool = False,
    separation: bool = False,
    pca_components: Optional[int] = None,
):
    """Run the full pipeline."""
    import anndata as ad
    import scanpy as sc
    
    model_name = get_model_name(checkpoint)
    hydra_overrides = get_hydra_overrides(checkpoint)
    
    panels_map = {
        "atlas": atlas_panels,
        "atlas2": atlas2_panels,
        "atlas_train": atlas_train_panels,
        "atlas2_train": atlas2_train_panels,
    }
    
    # Handle 'all' and 'all_all' options
    include_normal = {}
    for ds_name in ["atlas", "atlas2", "atlas_train", "atlas2_train"]:
        panels = panels_map.get(ds_name)
        include_normal[ds_name] = False
        if panels and "all_all" in panels:
            panels_map[ds_name] = list(PANELS.keys())
            include_normal[ds_name] = False  # all_all = only panels, no normal
        elif panels and "all" in panels:
            panels_map[ds_name] = ["zhuang", "zeng", "isd"]
            include_normal[ds_name] = True
    
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
    
    # Validate panels
    for ds_name, panels in panels_map.items():
        if panels and ds_name not in datasets:
            print(f"⚠️ Warning: --{ds_name.replace('_', '-')}-panel specified but '{ds_name}' not in datasets.")
            panels_map[ds_name] = None
            include_normal[ds_name] = False
    
    atlas_panels = panels_map["atlas"]
    atlas2_panels = panels_map["atlas2"]
    atlas_train_panels = panels_map["atlas_train"]
    atlas2_train_panels = panels_map["atlas2_train"]
    
    for ds_name, panels in panels_map.items():
        if panels:
            for panel in panels:
                if panel not in PANELS:
                    raise ValueError(f"Unknown panel: {panel}. Available: {list(PANELS.keys())}")
    
    def get_panels_for_ds(ds):
        return panels_map.get(ds)
    
    def get_include_normal(ds):
        return include_normal.get(ds, False)
    
    # Step 1: Generate embeddings
    print("=" * 60)
    print("STEP 1: EMBEDDINGS")
    print("=" * 60)
    
    for ds in datasets:
        ds_panels = get_panels_for_ds(ds)
        
        if is_atlas_dataset(ds) and ds_panels:
            for panel in ds_panels:
                emb_path = get_embedding_path(model_name, ds, panel=panel)
                
                if embeddings_exist(emb_path) and not force_regenerate:
                    print(f"✅ {ds}_{panel}: embeddings exist at {emb_path}")
                else:
                    filtered_atlas_path = get_filtered_atlas_path(ds, panel)
                    os.makedirs(os.path.dirname(filtered_atlas_path), exist_ok=True)
                    filter_adata_by_panel(
                        adata_path=CONFIG[ds]["adata_path"],
                        panel_path=PANELS[panel],
                        output_path=filtered_atlas_path,
                    )
                    generate_embeddings(
                        checkpoint=checkpoint,
                        adata_path=filtered_atlas_path,
                        output_path=emb_path,
                        hydra_overrides=hydra_overrides,
                        batch_size=batch_size,
                        subsample=subsample,
                    )
            
            if get_include_normal(ds):
                emb_path = get_embedding_path(model_name, ds)
                if embeddings_exist(emb_path) and not force_regenerate:
                    print(f"✅ {ds}: embeddings exist at {emb_path}")
                else:
                    generate_embeddings(
                        checkpoint=checkpoint,
                        adata_path=CONFIG[ds]["adata_path"],
                        output_path=emb_path,
                        hydra_overrides=hydra_overrides,
                        batch_size=batch_size,
                        subsample=subsample,
                    )
        else:
            emb_path = get_embedding_path(model_name, ds)
            if embeddings_exist(emb_path) and not force_regenerate:
                print(f"✅ {ds}: embeddings exist at {emb_path}")
            else:
                generate_embeddings(
                    checkpoint=checkpoint,
                    adata_path=CONFIG[ds]["adata_path"],
                    output_path=emb_path,
                    hydra_overrides=hydra_overrides,
                    batch_size=batch_size,
                    subsample=subsample,
                )
    
    # Step 2: Load datasets
    print("\n" + "=" * 60)
    print("STEP 2: LOAD DATASETS")
    print("=" * 60)
    
    adatas = {}
    for ds in datasets:
        ds_panels = get_panels_for_ds(ds)
        
        if is_atlas_dataset(ds) and ds_panels:
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
            emb_path = get_embedding_path(model_name, ds)
            adatas[ds] = load_dataset(
                dataset_name=ds,
                embedding_path=emb_path,
                subsample=CONFIG[ds].get("subsample"),
            )
        else:
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
        combined = ad.concat(adatas, join="outer", label="dataset", index_unique=None)
        combined.obs_names_make_unique()
    
    print(f"Combined: {combined.n_obs:,} cells | {dict(combined.obs['dataset'].value_counts())}")
    
    # Step 4: Compute neighbors
    mean_name = "mean" if mean else "cls"
    if compute_umap or clisi_keys or ilisi_keys:
        print("\n" + "=" * 60)
        print("STEP 4: NEIGHBORS + UMAP")
        print("=" * 60)
        
        use_rep = "concept_mean_embedding" if mean else "concept_cls_embedding"
        
        # Apply PCA if requested
        if pca_components:
            print(f"🔄 Applying PCA with {pca_components} components...")
            from sklearn.decomposition import PCA
            
            embeddings = combined.obsm[use_rep]
            pca = PCA(n_components=pca_components, random_state=0)
            combined.obsm["X_pca"] = pca.fit_transform(embeddings)
            
            explained_var = pca.explained_variance_ratio_.sum() * 100
            print(f"   PCA: {embeddings.shape[1]} dims -> {pca_components} dims ({explained_var:.1f}% variance)")
            
            use_rep = "X_pca"
            mean_name = f"{mean_name}_pca{pca_components}"
        
        sc.pp.neighbors(combined, use_rep=use_rep, n_neighbors=30)
    
    # Step 5: UMAP
    datasets_str = "_".join(datasets)  # default
    if compute_umap:
        sc.tl.umap(combined, min_dist=0.3, random_state=0)
        
        umap_dir = os.path.join(FIGURES_PATH, model_name)
        
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
        
        plot_umap(
            combined,
            color="dataset",
            title=f"UMAP by Dataset ({model_name})",
            out_png=os.path.join(umap_dir, f"{datasets_str}_dataset_{mean_name}.png"),
        )
        
        non_atlas_datasets = [ds for ds in datasets if not is_atlas_dataset(ds)]
        has_cell_type = any(CONFIG[ds]["has_cell_type"] for ds in non_atlas_datasets) if non_atlas_datasets else False
        if has_cell_type and "cell_type_mmc_raw" in combined.obs.columns:
            plot_umap(
                combined,
                color="cell_type_mmc_raw",
                palette=CELL_TYPE_PALETTE,
                title=f"UMAP by Cell Type ({model_name})",
                out_png=os.path.join(umap_dir, f"{datasets_str}_celltype_{mean_name}.png"),
            )
        
        # Separation: create separate UMAP for each cell type
        if separation and has_cell_type and "cell_type_mmc_raw" in combined.obs.columns:
            print("\n" + "=" * 60)
            print("STEP 4b: SEPARATION UMAPs (per cell type)")
            print("=" * 60)
            
            separation_dir = os.path.join(umap_dir, f"{datasets_str}_separation_{mean_name}")
            os.makedirs(separation_dir, exist_ok=True)
            
            cell_types = combined.obs["cell_type_mmc_raw"].dropna().unique()
            print(f"📊 Creating {len(cell_types)} separate UMAPs for each cell type...")
            
            for ct in sorted(cell_types):
                # Filter to this cell type
                mask = combined.obs["cell_type_mmc_raw"] == ct
                n_cells = mask.sum()
                
                if n_cells < 10:
                    print(f"   ⚠️ Skipping {ct}: only {n_cells} cells")
                    continue
                
                ct_adata = combined[mask].copy()
                
                # Plot colored by dataset
                ct_safe = ct.replace("/", "-").replace(" ", "_")
                plot_umap(
                    ct_adata,
                    color="dataset",
                    title=f"{ct} ({n_cells:,} cells)",
                    out_png=os.path.join(separation_dir, f"{ct_safe}_dataset.png"),
                )
                print(f"   ✅ {ct}: {n_cells:,} cells")
            
            print(f"\n💾 Separation UMAPs saved to: {separation_dir}")
    
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
    
    # Save results
    results_dir = os.path.join(FIGURES_PATH, model_name)
    os.makedirs(results_dir, exist_ok=True)
    results_file = os.path.join(results_dir, f"{datasets_str}_results_{mean_name}.txt")
    
    with open(results_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("EVALUATION RESULTS\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Model: {model_name}\n")
        f.write(f"Checkpoint: {checkpoint}\n")
        f.write(f"Datasets: {', '.join(datasets)}\n")
        if subsample:
            f.write(f"Subsample (embeddings): {subsample:,}\n")
        f.write("\n")
        for ds_name, panels in panels_map.items():
            if panels:
                f.write(f"{ds_name} panels: {', '.join(panels)}\n")
        f.write("\n" + "-" * 60 + "\n")
        f.write("METRICS\n")
        f.write("-" * 60 + "\n\n")
        for k, v in results.items():
            if k.startswith("cLISI") or k.startswith("iLISI"):
                f.write(f"{k}: {v:.4f}\n")
        f.write("\n" + "=" * 60 + "\n")
    
    print(f"\n💾 Results saved to: {results_file}")
    
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
    
    parser.add_argument("--checkpoint", "-c", required=True, help="Path to model checkpoint")
    parser.add_argument("--datasets", "-d", nargs="+", required=True,
                        choices=["isd", "zeng", "zhuang", "atlas", "atlas2", "atlas_train", "atlas2_train"],
                        help="Datasets to process")
    parser.add_argument("--umap", action="store_true", help="Compute and plot UMAP")
    parser.add_argument("--clisi", nargs="*", choices=["dataset", "cell_type"], help="Compute cLISI")
    parser.add_argument("--ilisi", nargs="*", choices=["dataset", "cell_type"], help="Compute iLISI")
    parser.add_argument("--batch-size", type=int, default=64, help="Batch size for embedding generation")
    parser.add_argument("--force", action="store_true", help="Force regeneration of embeddings")
    parser.add_argument("--atlas-panel", nargs="+",
                        help=f"Filter atlas to gene panel(s). 'all'=zhuang+zeng+isd+normal, 'all_all'=all panels in dir. Available: {list(PANELS.keys())}")
    parser.add_argument("--atlas2-panel", nargs="+",
                        help=f"Filter atlas2 to gene panel(s). Available: {list(PANELS.keys())}")
    parser.add_argument("--atlas-train-panel", nargs="+",
                        help=f"Filter atlas_train to gene panel(s). Available: {list(PANELS.keys())}")
    parser.add_argument("--atlas2-train-panel", nargs="+",
                        help=f"Filter atlas2_train to gene panel(s). Available: {list(PANELS.keys())}")
    parser.add_argument("--subsample", "-s", type=int, default=None,
                        help="Number of cells to subsample for embedding generation")
    parser.add_argument("--mean", action="store_true", help="Use mean embedding instead of CLS")
    parser.add_argument("--separation", action="store_true", 
                        help="Create separate UMAPs for each cell type, colored by dataset")
    parser.add_argument("--pca", type=int, default=None,
                        help="Apply PCA before UMAP (e.g., --pca 32 for 32 components)")
    
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
        subsample=args.subsample,
        mean=args.mean,
        separation=args.separation,
        pca_components=args.pca,
    )


if __name__ == "__main__":
    main()