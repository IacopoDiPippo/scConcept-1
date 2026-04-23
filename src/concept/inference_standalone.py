#!/usr/bin/env python3
"""
Standalone inference script for scConcept embeddings.
No dependencies on scConcept codebase, Hydra, or external configs.

Usage:
    python inference_standalone.py -c checkpoint.ckpt -g gene_mapping.pkl -i input.h5ad -o output.h5ad
    
Requirements:
    - torch
    - lightning
    - anndata
    - pandas
    - numpy
    - tqdm
"""

import argparse
import math
from pathlib import Path
from typing import Optional
from functools import partial

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
import anndata as ad


# ============================================================
# MODEL COMPONENTS (self-contained)
# ============================================================

class GeneEncoder(nn.Module):
    def __init__(self, n_genes: int, emb_dim: int, padding_idx: Optional[int] = None):
        super().__init__()
        self.embedding = nn.Embedding(n_genes, emb_dim, padding_idx=padding_idx)
        self.enc_norm = nn.LayerNorm(emb_dim)

    def forward(self, x: Tensor) -> Tensor:
        x = self.embedding(x)
        x = self.enc_norm(x)
        return x


class ContinuousValueEncoder(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.0):
        super().__init__()
        self.linear1 = nn.Linear(1, d_model)
        self.activation = nn.ReLU()
        self.linear2 = nn.Linear(d_model, d_model)
        self.norm = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, x: Tensor) -> Tensor:
        x = x.float().unsqueeze(-1)
        x = self.activation(self.linear1(x))
        x = self.linear2(x)
        x = self.norm(x)
        return self.dropout(x)


class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, max_len: int = 5000):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        seq_len = x.size(1)
        x = x + self.pe[:, :seq_len, :].to(x.device)
        return x


class StandaloneTransformerEncoder(nn.Module):
    """Simple transformer encoder without flash attention dependency."""
    
    def __init__(self, dim_model: int, num_head: int, dim_hid: int, nlayers: int, dropout: float = 0.1):
        super().__init__()
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim_model,
            nhead=num_head,
            dim_feedforward=dim_hid,
            dropout=dropout,
            batch_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=nlayers)
    
    def forward(self, x: Tensor, key_padding_mask: Tensor = None) -> Tensor:
        return self.encoder(x, src_key_padding_mask=key_padding_mask)


class StandaloneModel(nn.Module):
    """Standalone scConcept model for inference only."""
    
    CLS_VALUE = -2
    
    def __init__(
        self,
        dim_model: int,
        num_head: int,
        dim_hid: int,
        nlayers: int,
        vocab_size: int,
        pad_token_id: int,
        cls_token_id: int,
        dropout: float = 0.1,
        pe_max_len: int = 5000,
        projection_dim: Optional[int] = None,
        input_encoding: str = "rank_encoding",
        mask_padding: bool = True,
    ):
        super().__init__()
        
        self.dim_model = dim_model
        self.pad_token_id = pad_token_id
        self.cls_token_id = cls_token_id
        self.projection_dim = projection_dim
        self.input_encoding = input_encoding
        self.mask_padding = mask_padding
        
        # Encoders
        self.gene_token_encoder = GeneEncoder(vocab_size, dim_model, padding_idx=None)
        self.value_encoder = ContinuousValueEncoder(dim_model, dropout=0.0)
        self.positional_encoder = PositionalEncoding(dim_model, max_len=pe_max_len)
        
        # Transformer
        self.transformer_encoder = StandaloneTransformerEncoder(
            dim_model, num_head, dim_hid, nlayers, dropout
        )
        
        # Optional projection
        if projection_dim:
            self.projection = nn.Linear(dim_model, projection_dim, bias=False)
        else:
            self.projection = None
    
    def add_cls_token(self, tokens: Tensor, values: Tensor):
        batch_size = tokens.shape[0]
        
        cls_tokens = torch.full(
            (batch_size, 1), self.cls_token_id,
            dtype=tokens.dtype, device=tokens.device
        )
        cls_values = torch.full(
            (batch_size, 1), self.CLS_VALUE,
            dtype=values.dtype, device=values.device
        )
        
        tokens = torch.cat([cls_tokens, tokens], dim=1)
        values = torch.cat([cls_values, values], dim=1)
        
        return tokens, values
    
    def forward(self, tokens: Tensor, values: Tensor):
        # Add CLS token
        tokens, values = self.add_cls_token(tokens, values)
        
        # Padding mask
        if self.mask_padding:
            padding_mask = (tokens == self.pad_token_id)
        else:
            padding_mask = torch.zeros_like(tokens, dtype=torch.bool)
        
        # Encode genes
        gene_embs = self.gene_token_encoder(tokens)
        
        # Add positional or value encoding
        if self.input_encoding == 'rank_encoding':
            total_embs = self.positional_encoder(gene_embs)
        else:
            value_embs = self.value_encoder(values)
            total_embs = gene_embs + value_embs
        
        # Transform
        embs = self.transformer_encoder(total_embs, key_padding_mask=padding_mask)
        
        # CLS embedding
        cls_emb = embs[:, 0, :]
        
        # Mean embedding (excluding padding)
        mean_embs = []
        for i in range(len(embs)):
            mask = ~padding_mask[i]
            mean_embs.append(embs[i, mask].mean(dim=0))
        mean_emb = torch.stack(mean_embs, dim=0)
        
        # Optional projection
        if self.projection is not None:
            cls_emb = self.projection(cls_emb)
        
        return cls_emb, mean_emb


# ============================================================
# TOKENIZER
# ============================================================

class SimpleTokenizer:
    """Simple tokenizer for gene IDs."""
    
    def __init__(self, gene_mapping: dict):
        self.gene_to_token = gene_mapping
        self.pad_token = gene_mapping.get('<pad>', 0)
        self.cls_token = gene_mapping.get('<cls>', 1)
        self.vocab_size = len(gene_mapping)
    
    def tokenize(self, gene_ids: list) -> list:
        """Convert gene IDs to token IDs."""
        return [self.gene_to_token.get(g, self.pad_token) for g in gene_ids]


# ============================================================
# DATASET & COLLATE
# ============================================================

class SimpleDataset(Dataset):
    """Simple dataset for h5ad files."""
    
    def __init__(self, adata: ad.AnnData, tokenizer: SimpleTokenizer, gene_id_column: str = None):
        self.adata = adata
        self.tokenizer = tokenizer
        
        # Get gene IDs
        if gene_id_column and gene_id_column in adata.var.columns:
            self.gene_ids = adata.var[gene_id_column].tolist()
        else:
            self.gene_ids = adata.var_names.tolist()
        
        # Tokenize genes once
        self.gene_tokens = self.tokenizer.tokenize(self.gene_ids)
    
    def __len__(self):
        return self.adata.n_obs
    
    def __getitem__(self, idx):
        # Get expression values
        x = self.adata.X[idx]
        if hasattr(x, 'toarray'):
            x = x.toarray().flatten()
        else:
            x = np.array(x).flatten()
        
        return {
            'tokens': np.array(self.gene_tokens),
            'values': x.astype(np.float32),
        }


def collate_fn(batch, pad_token: int, max_tokens: int = 2048, gene_sampling_strategy: str = 'top-nonzero'):
    """Collate function with gene sampling."""
    
    tokens_list = []
    values_list = []
    
    for item in batch:
        tokens = item['tokens']
        values = item['values']
        
        # Get non-zero indices
        nonzero_mask = values > 0
        nonzero_indices = np.where(nonzero_mask)[0]
        
        if gene_sampling_strategy == 'top-nonzero':
            # Sort by value, take top max_tokens
            if len(nonzero_indices) > max_tokens:
                sorted_indices = nonzero_indices[np.argsort(values[nonzero_indices])[::-1]]
                selected = sorted_indices[:max_tokens]
            else:
                selected = nonzero_indices
        else:
            # Random sampling
            if len(nonzero_indices) > max_tokens:
                selected = np.random.choice(nonzero_indices, max_tokens, replace=False)
            else:
                selected = nonzero_indices
        
        # Sort by value (descending) for rank encoding
        selected = selected[np.argsort(values[selected])[::-1]]
        
        tokens_list.append(tokens[selected])
        values_list.append(values[selected])
    
    # Pad to max length in batch
    max_len = max(len(t) for t in tokens_list)
    
    padded_tokens = np.full((len(batch), max_len), pad_token, dtype=np.int64)
    padded_values = np.zeros((len(batch), max_len), dtype=np.float32)
    
    for i, (t, v) in enumerate(zip(tokens_list, values_list)):
        padded_tokens[i, :len(t)] = t
        padded_values[i, :len(v)] = v
    
    return {
        'tokens': torch.tensor(padded_tokens),
        'values': torch.tensor(padded_values),
    }


# ============================================================
# WEIGHT LOADING
# ============================================================

def load_model_from_checkpoint(checkpoint_path: str, gene_mapping: dict, device: torch.device):
    """Load model from Lightning checkpoint."""
    
    checkpoint = torch.load(checkpoint_path, map_location='cpu')
    
    # Get hyperparameters from checkpoint
    hparams = checkpoint.get('hyper_parameters', {})
    config = hparams.get('config', {})
    
    # Extract model config
    model_args = {
        'dim_model': config.get('dim_model', 512),
        'num_head': config.get('num_head', 8),
        'dim_hid': config.get('dim_hid', 1024),
        'nlayers': config.get('nlayers', 8),
        'dropout': config.get('dropout', 0.1),
        'pe_max_len': config.get('pe_max_len', 5000),
        'projection_dim': config.get('projection_dim', None),
        'input_encoding': config.get('input_encoding', 'rank_encoding'),
        'mask_padding': config.get('mask_padding', True),
        'vocab_size': len(gene_mapping),
        'pad_token_id': gene_mapping['<pad>'],
        'cls_token_id': gene_mapping['<cls>'],
    }
    
    print(f"📋 Model config:")
    for k, v in model_args.items():
        print(f"   {k}: {v}")
    
    # Create model
    model = StandaloneModel(**model_args)
    
    # Load state dict (need to map keys)
    state_dict = checkpoint['state_dict']
    
    # Map Lightning module keys to our model
    new_state_dict = {}
    for key, value in state_dict.items():
        # Remove 'model.' prefix if present
        new_key = key
        
        # Map transformer encoder keys
        if 'transformer_encoder.layers' in key:
            # Original: transformer_encoder.layers.0.xxx
            # Target: transformer_encoder.encoder.layers.0.xxx
            new_key = key.replace('transformer_encoder.layers', 'transformer_encoder.encoder.layers')
        
        new_state_dict[new_key] = value
    
    # Load with strict=False to handle any missing/extra keys
    missing, unexpected = model.load_state_dict(new_state_dict, strict=False)
    
    if missing:
        print(f"⚠️ Missing keys: {missing}")
    if unexpected:
        print(f"⚠️ Unexpected keys (ignored): {len(unexpected)} keys")
    
    model = model.to(device)
    model.eval()
    
    return model, model_args


# ============================================================
# MAIN INFERENCE FUNCTION
# ============================================================

def extract_embeddings(
    checkpoint_path: str,
    gene_mapping_path: str,
    input_path: str,
    output_path: str = None,
    batch_size: int = 64,
    max_tokens: int = 2048,
    gene_id_column: str = None,
    subsample: int = None,
):
    """Extract embeddings and add to h5ad file."""
    
    if output_path is None:
        output_path = input_path
    
    print(f"\n{'='*60}")
    print(f"🚀 Standalone Embedding Extraction")
    print(f"   Checkpoint: {checkpoint_path}")
    print(f"   Gene mapping: {gene_mapping_path}")
    print(f"   Input: {input_path}")
    print(f"   Output: {output_path}")
    print(f"{'='*60}\n")
    
    # Device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"📱 Device: {device}")
    
    # Load gene mapping
    print(f"📂 Loading gene mapping...")
    gene_mapping = pd.read_pickle(gene_mapping_path).to_dict()
    print(f"   Vocab size: {len(gene_mapping)}")
    
    # Load model
    print(f"🔧 Loading model...")
    model, model_args = load_model_from_checkpoint(checkpoint_path, gene_mapping, device)
    
    # Load adata
    print(f"📂 Loading {input_path}...")
    adata = ad.read_h5ad(input_path)
    print(f"   Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    
    # Subsample if requested
    if subsample and adata.n_obs > subsample:
        print(f"   Subsampling to {subsample:,} cells...")
        import scanpy as sc
        sc.pp.subsample(adata, n_obs=subsample, random_state=0)
    
    # Create tokenizer and dataset
    tokenizer = SimpleTokenizer(gene_mapping)
    dataset = SimpleDataset(adata, tokenizer, gene_id_column)
    
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=lambda b: collate_fn(b, tokenizer.pad_token, max_tokens),
        num_workers=0,
    )
    
    print(f"🔄 Extracting embeddings...")
    
    all_cls_embs = []
    all_mean_embs = []
    
    with torch.no_grad():
        with torch.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=device.type=='cuda'):
            for batch in tqdm(dataloader, desc="Generating embeddings"):
                tokens = batch['tokens'].to(device)
                values = batch['values'].to(device)
                
                cls_emb, mean_emb = model(tokens, values)
                
                all_cls_embs.append(cls_emb.float().cpu())
                all_mean_embs.append(mean_emb.float().cpu())
    
    # Concatenate
    cls_embeddings = torch.cat(all_cls_embs, dim=0).numpy()
    mean_embeddings = torch.cat(all_mean_embs, dim=0).numpy()
    
    print(f"   ✅ CLS embeddings: {cls_embeddings.shape}")
    print(f"   ✅ Mean embeddings: {mean_embeddings.shape}")
    
    # Add to adata
    adata.obsm["concept_cls_embedding"] = cls_embeddings
    adata.obsm["concept_mean_embedding"] = mean_embeddings
    
    # Save
    print(f"💾 Saving to {output_path}...")
    adata.write_h5ad(output_path)
    
    print(f"\n✅ Done!")
    
    return adata


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Standalone scConcept embedding extraction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument("--checkpoint", "-c", required=True, help="Path to model checkpoint (.ckpt)")
    parser.add_argument("--gene-mapping", "-g", required=True, help="Path to gene mapping pickle (.pkl)")
    parser.add_argument("--input", "-i", required=True, help="Path to input h5ad file")
    parser.add_argument("--output", "-o", default=None, help="Path to output h5ad file (default: overwrite input)")
    parser.add_argument("--batch-size", "-b", type=int, default=64, help="Batch size")
    parser.add_argument("--max-tokens", "-m", type=int, default=2048, help="Max tokens per cell")
    parser.add_argument("--gene-id-column", default=None, help="Column in adata.var for gene IDs (default: var_names)")
    parser.add_argument("--subsample", "-s", type=int, default=None, help="Subsample to N cells")
    
    args = parser.parse_args()
    
    extract_embeddings(
        checkpoint_path=args.checkpoint,
        gene_mapping_path=args.gene_mapping,
        input_path=args.input,
        output_path=args.output,
        batch_size=args.batch_size,
        max_tokens=args.max_tokens,
        gene_id_column=args.gene_id_column,
        subsample=args.subsample,
    )


if __name__ == "__main__":
    main()