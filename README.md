# scConcept

<!-- [![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/theislab/scConcept/test.yaml?branch=main
[badge-docs]: https://img.shields.io/readthedocs/scConcept -->

This repository contains the python package to train and use scConcept (Single-cell contrastive cell pre-training) method for single-cell transcriptomics.

<!-- ## Getting started

Please refer to the [documentation][],
in particular, the [API documentation][]. -->

## Installation

You need to have Python 3.10 or newer installed on your system.
If you don't have Python installed, we recommend installing [uv][].

scConcept also uses [Flash Attention][] which requires CUDA

<!-- There are several alternative options to install scConcept: -->

<!--
1) Install the latest release of `scConcept` from [PyPI][]:

```bash
pip install scConcept
```
-->

1. Install the latest development version:

```bash
pip install git+https://github.com/theislab/scConcept.git@main
```

2. [Flash Attention][] (required) - CUDA is required for installing flash-attn:

```bash
pip install flash-attn==2.7.* --no-build-isolation
```

3. Install lamin-dataloader from GitHub (required):

```bash
pip install git+https://github.com/theislab/lamin_dataloader.git
```

# scConcept_MOUSE

Mouse brain adaptation of [scConcept](https://github.com/theislab/scConcept) — a contrastive learning framework for technology-agnostic single-cell embeddings.

## Overview

This project adapts scConcept to learn cell embeddings from **mouse brain** single-cell and spatial transcriptomics data. The model is pretrained on dissociated scRNA-seq data (ABC Atlas) using a contrastive objective, producing embeddings that generalize across different spatial gene panels and technologies. These embeddings serve as input features for a downstream **point transformer** model for spatial cell analysis.

### What's in this repo

- **Mouse adaptation** of scConcept's contrastive pretraining pipeline
- **Weighted gene sampling** based on differential expression (DE), PCA loadings, and highly variable genes (HVG)
- **Finetuning strategies** on spatial MERFISH data (Zeng, Zhuang, ISD)
- **Curated mouse brain gene panels** from MERFISH/MERSCOPE codebooks
- **Evaluation pipeline**: embedding generation, UMAP visualization, cLISI/iLISI metrics
- **Cell type prediction** benchmarks on embeddings

## Project Structure
```
scConcept_MOUSE/
├── src/
│   └── concept/
│       ├── train.py              # Main training script
│       ├── evaluate.py           # Evaluation pipeline (embeddings, UMAP, metrics)
│       ├── get_embs.py           # Embedding generation module
│       └── ...
├── configs/
│   ├── model/
│   │   └── ContrastiveModel.yaml
│   ├── datamodule/
│   │   └── DataModuleBasic.yaml
│   └── split/
│       ├── split_mouse.yaml
│       └── split_mouse_2.yaml
├── Panels/
│   ├── panels_with_full_names/       # Panel CSVs (original names)
│   ├── panels_with_shortened_names/  # Panel CSVs (short names)
│   ├── various_panels/               # token_mapping, sbatch examples, overrides
│   └── panel_with_genes_scores.csv   # Per-gene weights for weighted sampling
├── notebooks/
│   └── cell_type_pred.ipynb      # Cell type classification on embeddings
├── handover.pdf                  # Full project documentation
└── README.md
```

## Setup

### Environment (JURECA)

1. Create a virtual environment following the [JSC guide](https://jupyterjsc.pages.jsc.fz-juelich.de/docs/jupyterjsc/users/jupyterlab/4.3/kernels_hpc_venv/).
2. Install dependencies from the [original scConcept repo](https://github.com/theislab/scConcept).
3. Log in to Weights & Biases: `wandb login`

There is an activation script (`.sh`) in the repo for convenience.

### Data Location (JURECA)
```
/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/concept_embeddings/
├── train_val_data/       # .h5ad files for training and validation
├── atlas_filtered/       # Atlas filtered to specific spatial panels
├── models/               # Generated embeddings (by model and dataset)
├── important/finetune/   # Finetuned model checkpoints
├── hvg_isd_normed.h5ad   # ISD dataset (normalized)
└── ISD-1.h5ad            # ISD dataset (raw)
```

## Training

Submit jobs via `sbatch` (example scripts in `Panels/various_panels/`). Configuration is managed through Hydra YAML files — override on the command line:
```bash
python -m concept.train \
    model=ContrastiveModel \
    datamodule=DataModuleBasic \
    split=split_mouse \
    wandb.run_name=my_experiment \
    datamodule.probabilistic_panel_sampling=True
```

> **Note:** GPU nodes on JURECA have no internet. W&B logs offline. After training, sync from a login node:
> ```bash
> wandb sync offline-run*
> ```

## Evaluation

Run from the `src/` directory:
```bash
# Full evaluation with UMAP and metrics
python -m concept.evaluate -c /path/to/checkpoint.ckpt \
    -d isd zeng zhuang --umap --clisi cell_type --ilisi dataset

# Atlas with spatial panel views
python -m concept.evaluate -c /path/to/checkpoint.ckpt \
    -d atlas --atlas-panel all --umap --ilisi dataset

# PCA 32 (point transformer input dimensionality)
python -m concept.evaluate -c /path/to/checkpoint.ckpt \
    -d isd zeng zhuang --umap --pca 32

# Per-cell-type integration analysis
python -m concept.evaluate -c /path/to/checkpoint.ckpt \
    -d isd zeng zhuang --umap --separation
```

See `handover.pdf` for the full flag reference and all usage examples.

## Trained Models

### Baseline (pretrained on transcriptomics only)

| Atlas | Size | Sampling | W&B Run ID |
|-------|------|----------|------------|
| Atlas 1 | Big | Uniform | `z6csa92l` |
| Atlas 1 | Big | Weighted | `uv28i4o3` |
| Atlas 1 | Small | Uniform | `ko0hq20i` |
| Atlas 1 | Small | Weighted | `kdp2po4o` |
| Atlas 2 | Big | Uniform | `b5wnm0b2` |
| Atlas 2 | Big | Weighted | `6iz4xrqm` |
| Atlas 2 | Small | Uniform | `dfa4zrfz` |
| Atlas 2 | Small | Weighted | `hgld4ax1` |

### Finetuned (from small weighted model)

| Checkpoint | Description |
|------------|-------------|
| `3ektavyt` | Panel exposure only (no spatial data) |
| `3de8fprx` | Spatial finetuning: Zeng + Zhuang |
| `2cf8i4ws` | Spatial finetuning: Zeng + Zhuang + ISD |
| `2i8eeuo4` | Spatial finetuning: Zeng + Zhuang + ISD (second run with also extra panel ( first panel always determined, second with 25%) |
| `n1emty30` | Other finetuning experiments |
| `n1g2uu4t` | Other finetuning experiments |
| `vrxkd74x` | Other finetuning experiments |

Checkpoints stored at: `concept_embeddings/important/finetune/`

## Gene Panels

Curated from MERFISH/MERSCOPE codebooks (mouse brain only). All panels normalized to Ensembl gene IDs (ENSMUSG). Panels excluded from training to avoid redundancy: `Allen_Zeng` (≈identical to `Zhuang_ABCA1`) and `MouseBrain_Final` (≈identical to `ISD`).

See `Panels/` directory and `handover.pdf` for full panel documentation and overlap analysis.

## Documentation

The complete project handover document is in `handover.pdf`. It covers:

- Scientific background (omics data, contrastive learning, scConcept architecture)
- All modifications made (mouse adaptation, weighted sampling, finetuning)
- Results with UMAPs and metrics for every experiment
- Full YAML configuration reference with explanation of every field
- Complete `evaluate.py` command reference
- Data and checkpoint locations

## References

- [scConcept (original)](https://github.com/theislab/scConcept) — Theislab, 2025
- [scConcept preprint](https://www.biorxiv.org/content/10.1101/2025.10.14.682419v1) — bioRxiv, 2025



## Citation

> Bahrami, M., Tejada-Lapuerta, A., Becker, S., Hashemi G, F.S. and Theis, F.J., 2025. scConcept: Contrastive pretraining for technology-agnostic single-cell representations beyond reconstruction. bioRxiv, pp.2025-10. doi: https://doi.org/10.1101/2025.10.14.682419

[uv]: https://github.com/astral-sh/uv
[Flash Attention]: https://github.com/Dao-AILab/flash-attention
[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/theislab/scConcept/issues