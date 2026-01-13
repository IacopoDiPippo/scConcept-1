import os
import sys
import torch
from omegaconf import DictConfig, OmegaConf
from concept import scConcept
import numpy as np
from pathlib import Path
import argparse
import anndata as ad
from hydra import compose, initialize


def get_embs(cfg: DictConfig, ckpt_path: str, adata_path: str, gene_id_column: str = None,
             batch_size: int = 32, max_tokens: int = None, gene_sampling_strategy: str = None):
    """
    Extract embeddings using scConcept API extract_embeddings method.
    
    Args:
        cfg: Configuration dictionary
        ckpt_path: Path to model checkpoint
        adata_path: Path to AnnData file
        batch_size: Batch size for dataloader
        max_tokens: Maximum number of tokens per cell
        gene_sampling_strategy: Gene sampling strategy
        gene_id_column: Column name in adata.var to use as gene IDs (default: None, uses index)
    """
    concept = scConcept()
    concept.load_config_and_model(
        config=cfg,
        model_path=ckpt_path,
        gene_mapping_path=cfg.PATH.gene_mapping_path,
        panels_dir=cfg.PATH.PANELS_PATH
    )
    
    print(f"Loading AnnData from {adata_path}...")
    adata = ad.read_h5ad(adata_path)
        
    print(f"Extracting embeddings with batch_size={batch_size}, max_tokens={max_tokens}, gene_sampling_strategy={gene_sampling_strategy}")
    result = concept.extract_embeddings(
        adata=adata,
        batch_size=batch_size,
        max_tokens=max_tokens,
        gene_sampling_strategy=gene_sampling_strategy,
        gene_id_column=gene_id_column
    )
    
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=str, required=True, help="Path to checkpoint file")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to the AnnData file (.h5ad)")
    parser.add_argument("--gene_id_column", type=str, default=None, help="Column name in adata.var to use as gene IDs")
    parser.add_argument("--output_emb_path", type=str, required=True, help="Path to save embeddings")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--max_tokens", type=int, default=None, help="Maximum tokens per cell")
    parser.add_argument("--gene_sampling_strategy", type=str, default=None, help="Gene sampling strategy")
    args, unknown = parser.parse_known_args()
    
    print(f"Loading config from Hydra...")
    print('overrides:', unknown)
    with initialize(version_base=None, config_path="./conf"):
        cfg = compose(config_name="config", overrides=unknown)
    
    gpu_info = torch.cuda.get_device_properties(0)
    print(f"GPU Type: {gpu_info.name}")
    
    output_emb_path = Path(args.output_emb_path)
    
    if (output_emb_path / 'cell_embs_cls.npy').exists():
        print(f"Embeddings already exist in {output_emb_path}")
        exit(0)
    
    print(f"Using checkpoint: {args.checkpoint}")
    
    result = get_embs(
        cfg=cfg, 
        ckpt_path=args.checkpoint, 
        adata_path=args.adata_path,
        gene_id_column=args.gene_id_column,
        batch_size=args.batch_size,
        max_tokens=args.max_tokens,
        gene_sampling_strategy=args.gene_sampling_strategy,
    )
    
    print(f"Saving embeddings to {output_emb_path}...")
    os.makedirs(output_emb_path, exist_ok=True)
    np.save(output_emb_path / 'cell_embs_cls.npy', result['cls_cell_emb'])
    np.save(output_emb_path / 'cell_embs_mean.npy', result['mean_cell_emb'])
    
    if 'context_sizes' in result:
        np.save(output_emb_path / 'context_sizes.npy', result['context_sizes'])
    
    print(f"Embeddings saved successfully to {output_emb_path}")



"""
python src/concept/get_embs.py \
  --checkpoint /p/project1/hai_fzj_bda/checkpoints/t9qa3400/steps/step=310000.ckpt \
  --adata_path /p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/my_dataset.h5ad \
  --output_emb_path /p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/embeddings/my_dataset \
  --batch_size 64
  """