#!/usr/bin/env python
"""
Cell Type Classification Script

Usage:
    python celltype_classification.py --dataset zeng --model big_weighted__uv28i4o3_last --label_key class
"""

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn.functional as F
from sklearn.model_selection import train_test_split
from ct_rep.celltype.celltype_pred import evaluate

# Dataset configurations
DATASET_CONFIG = {
    "zeng": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zeng.h5ad",
        "dataname": "zeng",
        "needs_split": True,
        "annotation_path": "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zeng/results/merfish/cell_type_annotation/adata_obs_annotated.csv",
        "cell_id_suffix": "-Zeng",
    },
    "zhuang": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zhuang-ABCA-1.h5ad",
        "dataname": "zhuang",
        "needs_split": True,
        "annotation_path": "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv",
        "cell_id_suffix": "-Zhuang-ABCA-1",
    },
    "isd": {
        "adata_path": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/hvg_isd_normed.h5ad",
        "dataname": "isd",
        "needs_split": True,
        "annotation_path": None,  # ISD has cell_type in adata already
        "cell_id_suffix": None,
    },
    "atlas": {
        "adata_path_train": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/train_val_data/abc_atlas_train.h5ad",
        "adata_path_test": "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/train_val_data/abc_atlas_val.h5ad",
        "dataname_train": "atlas_train",
        "dataname_test": "atlas",
        "needs_split": False,
        "annotation_path": None,
        "cell_id_suffix": None,
    },
}

# Model ID -> (size, weighting, atlas_version)
# split_mouse_2 = atlas2, others = atlas1
MODEL_INFO = {
    # Big models - split_mouse_2 (atlas2)
    "pzrmdo2r": ("big", "uniform", "atlas2"),
    "qisq633l": ("big", "weighted", "atlas2"),
    "b5wnm0b2": ("big", "uniform", "atlas2"),
    "6iz4xrqm": ("big", "weighted", "atlas2"),
    # Small models - split_mouse_2 (atlas2)
    "t2hxc3hd": ("small", "uniform", "atlas2"),
    "7yxo1joa": ("small", "weighted", "atlas2"),
    "dfa4zrfz": ("small", "uniform", "atlas2"),
    "hgld4ax1": ("small", "weighted", "atlas2"),
}

EMB_BASE_PATH = "/p/home/jusers/dipippo1/jureca/projects/scConcept-1/src/concept/concept_embeddings/models"


def parse_model_id(model_id):
    """Parse model_id to extract size, weighting, and atlas version."""
    # Extract the run ID from model_id (e.g., "big_weighted__uv28i4o3_last" -> "uv28i4o3")
    for run_id in MODEL_INFO:
        if run_id in model_id:
            return MODEL_INFO[run_id]
    
    # Fallback: try to parse from string
    model_lower = model_id.lower()
    size = "big" if "big" in model_lower else ("small" if "small" in model_lower else "unknown")
    weighting = "weighted" if "weighted" in model_lower else ("uniform" if "uniform" in model_lower else "unknown")
    atlas_version = "atlas1"  # default
    
    return size, weighting, atlas_version


def add_cell_type_mmc_raw(adata, dataset: str, dataset_cfg: dict):
    """
    Add cell_type_mmc_raw annotation to adata.
    Only works for zeng, zhuang, isd.
    """
    annotation_path = dataset_cfg.get("annotation_path")
    cell_id_suffix = dataset_cfg.get("cell_id_suffix")
    
    if dataset == "isd":
        # ISD has cell_type in adata already
        if "cell_type_mmc_raw_revised" in adata.obs.columns:
            adata.obs["cell_type_mmc_raw"] = adata.obs["cell_type_mmc_raw_revised"]
        elif "cell_type" in adata.obs.columns:
            adata.obs["cell_type_mmc_raw"] = adata.obs["cell_type"]
        else:
            raise ValueError("ISD: No cell_type column found")
        
        # Filter to annotated cells
        mask_valid = adata.obs["cell_type_mmc_raw"].notna()
        print(f"   Annotated: {mask_valid.sum()}/{adata.n_obs} cells")
        adata = adata[mask_valid].copy()
        
    elif annotation_path is not None:
        # Zeng and Zhuang - load from CSV
        ann = pd.read_csv(annotation_path)
        
        # Create cell_id with proper suffix
        adata.obs["cell_id"] = adata.obs.index.astype(str) + cell_id_suffix
        
        # For Zeng, also add suffix to annotation CSV (it doesn't have it)
        if dataset == "zeng":
            ann["cell_id"] = ann["cell_id"].astype(str) + cell_id_suffix
        
        # Merge annotations
        adata.obs = adata.obs.merge(
            ann[["cell_id", "cell_type_mmc_raw"]],
            on="cell_id",
            how="left"
        )
        adata.obs.index = adata.obs["cell_id"]
        
        # Filter to annotated cells only
        mask_valid = adata.obs["cell_type_mmc_raw"].notna()
        print(f"   Annotated: {mask_valid.sum()}/{adata.n_obs} cells")
        adata = adata[mask_valid].copy()
    else:
        raise ValueError(f"Dataset {dataset} does not support cell_type_mmc_raw")
    
    return adata


def load_embeddings(adata, model_id, dataname):
    """Load embeddings for a given model and align with adata."""
    emb_path = f"{EMB_BASE_PATH}/{model_id}/{dataname}/cell_embs_cls.npy"
    cell_ids_path = f"{EMB_BASE_PATH}/{model_id}/{dataname}/cell_ids.npy"
    
    if not os.path.exists(emb_path):
        raise FileNotFoundError(f"Embedding file not found: {emb_path}")
    
    emb = np.load(emb_path)
    
    if os.path.exists(cell_ids_path):
        cell_ids = np.load(cell_ids_path, allow_pickle=True).astype(str)
        common_cells = adata.obs_names.intersection(cell_ids)
        print(f"Cells in adata: {adata.n_obs}, Cells with embeddings: {len(cell_ids)}, Common: {len(common_cells)}")
        
        if len(common_cells) == 0:
            raise ValueError("No common cells found between adata and embeddings")
        
        adata = adata[common_cells].copy()
        emb_df = pd.DataFrame(emb, index=cell_ids)
        adata.obsm['X_emb'] = emb_df.loc[adata.obs_names].values
    else:
        adata.obsm['X_emb'] = emb
    
    return adata


def run_classification(
    dataset: str,
    model_id: str,
    label_key: str,
    classifier: str = "knn",
    min_count: int = 50,
    train_size: int = 200000,
    test_size: int = 60000,
    test_split: float = 0.2,
    results_path: str = "classify_results.csv",
):
    if dataset not in DATASET_CONFIG:
        raise ValueError(f"Unknown dataset: {dataset}. Choose from {list(DATASET_CONFIG.keys())}")
    
    # cell_type_mmc_raw only works for zeng, zhuang, isd
    if label_key == "cell_type_mmc_raw" and dataset not in ["zeng", "zhuang", "isd"]:
        raise ValueError(f"cell_type_mmc_raw only works for zeng, zhuang, isd. Got: {dataset}")
    
    dataset_cfg = DATASET_CONFIG[dataset]
    size, weighting, atlas_version = parse_model_id(model_id)
    
    print(f"\nDataset: {dataset} | Model: {model_id} | Label: {label_key}")
    print(f"Size: {size} | Weighting: {weighting} | Atlas: {atlas_version}\n")
    
    # Load data
    if dataset_cfg["needs_split"]:
        adata = sc.read_h5ad(dataset_cfg["adata_path"])
        adata = load_embeddings(adata, model_id, dataset_cfg["dataname"])
        
        # Add cell_type_mmc_raw if needed
        if label_key == "cell_type_mmc_raw":
            print("Loading cell_type_mmc_raw annotations...")
            adata = add_cell_type_mmc_raw(adata, dataset, dataset_cfg)
        
        # Filter low frequency
        freq = adata.obs[label_key].value_counts()
        adata = adata[adata.obs[label_key].isin(freq[freq >= min_count].index)].copy()
        
        # Split
        indices = np.arange(adata.n_obs)
        train_idx, test_idx = train_test_split(indices, test_size=test_split, random_state=42, stratify=adata.obs[label_key])
        adata_train = adata[train_idx].copy()
        adata_test = adata[test_idx].copy()
    else:
        adata_train = sc.read_h5ad(dataset_cfg["adata_path_train"])
        adata_test = sc.read_h5ad(dataset_cfg["adata_path_test"])
        
        # Load embeddings
        emb_path_train = f"{EMB_BASE_PATH}/{model_id}/{dataset_cfg['dataname_train']}/cell_embs_cls.npy"
        emb_path_test = f"{EMB_BASE_PATH}/{model_id}/{dataset_cfg['dataname_test']}/cell_embs_cls.npy"
        
        if not os.path.exists(emb_path_train) or not os.path.exists(emb_path_test):
            raise FileNotFoundError(f"Embedding files not found")
        
        adata_train.obsm['X_emb'] = np.load(emb_path_train)
        adata_test.obsm['X_emb'] = np.load(emb_path_test)
        
        # Filter low frequency
        freq = adata_train.obs[label_key].value_counts()
        adata_train = adata_train[adata_train.obs[label_key].isin(freq[freq >= min_count].index)].copy()
        adata_test = adata_test[adata_test.obs[label_key].isin(freq[freq >= min_count].index)].copy()
    
    # Subsample
    if train_size and adata_train.n_obs > train_size:
        sc.pp.subsample(adata_train, n_obs=train_size, random_state=42)
    if test_size and adata_test.n_obs > test_size:
        sc.pp.subsample(adata_test, n_obs=test_size, random_state=42)
    
    # Normalize (cosine)
    adata_train.obsm['X_emb'] = F.normalize(torch.tensor(adata_train.obsm['X_emb']), p=2, dim=1).numpy()
    adata_test.obsm['X_emb'] = F.normalize(torch.tensor(adata_test.obsm['X_emb']), p=2, dim=1).numpy()
    
    print(f"Train: {adata_train.n_obs} cells | Test: {adata_test.n_obs} cells")
    
    # Evaluate
    accuracy, f1_weighted, f1_macro, _ = evaluate(
        adata_train.obsm['X_emb'], adata_train.obs[label_key],
        adata_test.obsm['X_emb'], adata_test.obs[label_key],
        classifier, metric="cosine", pred_within_test_labels=False, lr=0.0005
    )
    
    print(f"Accuracy: {accuracy:.4f} | F1 Weighted: {f1_weighted:.4f} | F1 Macro: {f1_macro:.4f}")
    
    # Save result
    result = {
        'dataset': dataset,
        'model': model_id,
        'label_key': label_key,
        'size': size,
        'weighting': weighting,
        'atlas_version': atlas_version,
        'accuracy': accuracy,
        'f1_weighted': f1_weighted,
        'f1_macro': f1_macro,
    }
    
    columns = ['dataset', 'model', 'label_key', 'size', 'weighting', 'atlas_version', 'accuracy', 'f1_weighted', 'f1_macro']
    
    if os.path.exists(results_path):
        results_df = pd.read_csv(results_path)
    else:
        results_df = pd.DataFrame(columns=columns)
    
    # Append (no duplicate removal - keep all runs)
    results_df = pd.concat([results_df, pd.DataFrame([result])], ignore_index=True)
    results_df.to_csv(results_path, index=False)
    print(f"Saved to {results_path}")
    
    return result


def main():
    parser = argparse.ArgumentParser(description="Cell Type Classification")
    parser.add_argument("--dataset", type=str, required=True, choices=["zeng", "zhuang", "isd", "atlas"])
    parser.add_argument("--model", type=str, required=True, help="Model ID (e.g., big_weighted__uv28i4o3_last)")
    parser.add_argument("--label_key", type=str, required=True, help="Label column (e.g., class)")
    parser.add_argument("--classifier", type=str, default="knn")
    parser.add_argument("--min_count", type=int, default=50)
    parser.add_argument("--train_size", type=int, default=200000)
    parser.add_argument("--test_size", type=int, default=60000)
    parser.add_argument("--results_path", type=str, default="classify_results.csv")
    
    args = parser.parse_args()
    
    run_classification(
        dataset=args.dataset,
        model_id=args.model,
        label_key=args.label_key,
        classifier=args.classifier,
        min_count=args.min_count,
        train_size=args.train_size,
        test_size=args.test_size,
        results_path=args.results_path,
    )


if __name__ == "__main__":
    main()