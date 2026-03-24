#!/usr/bin/env python3
"""
Add scConcept embeddings to an h5ad file.

Usage:
    python -m concept.add_embeddings -c /path/to/checkpoint.ckpt -i input.h5ad -o output.h5ad
    
    # Overwrite input file
    python -m concept.add_embeddings -c /path/to/checkpoint.ckpt -i input.h5ad
    
    # With subsample
    python -m concept.add_embeddings -c /path/to/checkpoint.ckpt -i input.h5ad -o output.h5ad -s 50000
"""

import argparse
import os
from pathlib import Path
from typing import List

import numpy as np
import anndata as ad


def get_checkpoint_dir(checkpoint_path: str) -> Path:
    """Get checkpoint directory containing overrides.txt."""
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
    skip_prefixes = (
        'wandb.',
        '+initialize_2.',
        'datamodule.dataloader.',
        'datamodule.probabilistic_panel_sampling',
    )
    filtered = [o for o in overrides if not any(o.startswith(p) for p in skip_prefixes)]
    
    print(f"📋 Loaded overrides from {overrides_file}:")
    for o in filtered:
        print(f"   {o}")
    
    return filtered


def add_embeddings(
    checkpoint: str,
    input_path: str,
    output_path: str = None,
    batch_size: int = 64,
    subsample: int = None,
):
    """
    Add scConcept embeddings to h5ad file.
    
    Args:
        checkpoint: Path to model checkpoint
        input_path: Path to input h5ad file
        output_path: Path to output h5ad file (if None, overwrites input)
        batch_size: Batch size for inference
        subsample: Number of cells to subsample (None for all)
    """
    import scanpy as sc
    from omegaconf import DictConfig
    from hydra import compose, initialize_config_dir, initialize
    import torch
    
    from concept.api import scConcept
    
    if output_path is None:
        output_path = input_path
    
    print(f"\n{'='*60}")
    print(f"🚀 Adding embeddings to h5ad")
    print(f"   Checkpoint: {checkpoint}")
    print(f"   Input: {input_path}")
    print(f"   Output: {output_path}")
    if subsample:
        print(f"   Subsample: {subsample:,} cells")
    print(f"{'='*60}\n")
    
    # Load overrides
    hydra_overrides = get_hydra_overrides(checkpoint)
    
    # Initialize Hydra and load config
    print("📋 Loading config...")
    with initialize(version_base=None, config_path="./conf"):
        cfg = compose(config_name="config", overrides=hydra_overrides)
    
    # Load model
    print("🔧 Loading model...")
    concept = scConcept()
    concept.load_config_and_model(
        config=cfg,
        model_path=checkpoint,
        gene_mapping_path=cfg.PATH.gene_mapping_path,
        panels_dir=cfg.PATH.PANELS_PATH
    )
    
    # Load adata
    print(f"📂 Loading {input_path}...")
    adata = ad.read_h5ad(input_path)
    print(f"   Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    
    # Subsample if requested
    if subsample and adata.n_obs > subsample:
        print(f"   Subsampling to {subsample:,} cells...")
        sc.pp.subsample(adata, n_obs=subsample, random_state=0)
    
    # Extract embeddings
    print(f"🔄 Extracting embeddings (batch_size={batch_size})...")
    result = concept.extract_embeddings(
        adata=adata,
        batch_size=batch_size,
    )
    
    # Add embeddings to adata
    adata.obsm["concept_cls_embedding"] = result['cls_cell_emb']
    adata.obsm["concept_mean_embedding"] = result['mean_cell_emb']
    
    print(f"   ✅ Added concept_cls_embedding: {result['cls_cell_emb'].shape}")
    print(f"   ✅ Added concept_mean_embedding: {result['mean_cell_emb'].shape}")
    
    # Save
    print(f"💾 Saving to {output_path}...")
    adata.write_h5ad(output_path)
    
    print(f"\n✅ Done! Embeddings added to {output_path}")
    
    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Add scConcept embeddings to h5ad file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    
    parser.add_argument("--checkpoint", "-c", required=True, help="Path to model checkpoint")
    parser.add_argument("--input", "-i", required=True, help="Path to input h5ad file")
    parser.add_argument("--output", "-o", default=None, help="Path to output h5ad file (default: overwrite input)")
    parser.add_argument("--batch-size", "-b", type=int, default=64, help="Batch size for inference")
    parser.add_argument("--subsample", "-s", type=int, default=None, help="Subsample to N cells")
    
    args = parser.parse_args()
    
    add_embeddings(
        checkpoint=args.checkpoint,
        input_path=args.input,
        output_path=args.output,
        batch_size=args.batch_size,
        subsample=args.subsample,
    )


if __name__ == "__main__":
    main()