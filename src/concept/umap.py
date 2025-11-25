import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------
# Helpers
# -----------------------------------------
def neighbors_umap(adata, use_rep=None):
    print("PP Neighbors")
    if use_rep:
        sc.pp.neighbors(adata, use_rep=use_rep)
    else:
        sc.pp.neighbors(adata)

    print("UMAP")
    sc.tl.umap(adata, min_dist=0.3, random_state=0)


def plot(adata, color, out_png=None, title=None, palette=None):
    if title is None:
        title = f"UMAP colored by {color}"

    if palette:
        ax = sc.pl.umap(
            adata,
            color=color,
            title=title,
            frameon=False,
            palette=palette,
            show=False,
            return_fig=True,
        )
    else:
        ax = sc.pl.umap(
            adata,
            color=color,
            title=title,
            frameon=False,
            show=False,
            return_fig=True,
        )

    plt.show()

    if out_png is not None:
        print(f"Saving figure to: {out_png}")
        ax.savefig(out_png, dpi=300, bbox_inches="tight")

    plt.close(ax)


def sanity_check(adata, color=None):
    print("\n======================")
    print("🔍   ANNDATA CHECK")
    print("======================")
    print(f"Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print("obs columns:", list(adata.obs.columns))
    print("obsm keys:", list(adata.obsm.keys()))

    if color is not None and color in adata.obs.columns:
        missing = adata.obs[color].isna().sum()
        print(f"\nColumn '{color}':")
        print(f"  unique values: {adata.obs[color].nunique()}")
        print(f"  missing values: {missing:,}")
    else:
        if color is not None:
            print(f"\n⚠️ Column '{color}' not found in adata.obs")

    if "X_umap" in adata.obsm:
        print("UMAP present, shape:", adata.obsm["X_umap"].shape)
    else:
        print("No UMAP yet (X_umap missing).")
    print("======================\n")


# -----------------------------------------
# Fixed palette
# -----------------------------------------
cell_type_palette = {
    'ABCs': '#023fa5',
    'Astrocytes': '#7d87b9',
    'Astroependymal': '#bec1d4',
    'BAMs': '#d6bcc0',
    'Bergmann': '#bb7784',
    'Choroid-Plexus': '#8e063b',
    'ECs': '#4a6fe3',
    'Ependymal': '#8595e1',
    'Immune-Other': '#b5bbe3',
    'Microglia': '#e6afb9',
    'Neurons-Dopa': '#e07b91',
    'Neurons-Gaba': '#d33f6a',
    'Neurons-Glut': '#11c638',
    'Neurons-Glyc-Gaba': '#8dd593',
    'Neurons-Granule-Immature': '#c6dec7',
    'Neurons-Other': '#ead3c6',
    'OECs': '#f0b98d',
    'OPCs': '#ef9708',
    'Oligodendrocytes': '#0fcfc0',
    'Pericytes': '#9cded6',
    'SMCs': '#d5eae7',
    'Tanycytes': '#f3e1eb',
    'Undefined': '#f6c4e1',
    'VLMCs': '#f79cd4',
}


# -----------------------------------------
# 1) Load AnnData with concept embeddings
# -----------------------------------------
adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zeng_concept_embeddings.h5ad"
print(f"📂 Loading: {adata_path}")
adata = sc.read(adata_path)

print("🧪 Sanity check BEFORE anything:")
sanity_check(adata)


# -----------------------------------------
# 2) Add cell_id and subsample (like zhuang code)
# -----------------------------------------
adata.obs["cell_id"] = adata.obs.index

print("📉 Subsampling to 100,000 cells (if available)...")
sc.pp.subsample(adata, n_obs=100000, random_state=0)

print("🧪 Sanity check AFTER subsample (still before annotations):")
sanity_check(adata)


# -----------------------------------------
# 3) Attach annotations (Zhuang-style)
# -----------------------------------------
annotations_file = "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv"
print(f"📂 Loading annotations from: {annotations_file}")
ann = pd.read_csv(annotations_file)

adata.obs = adata.obs.merge(
    ann[["cell_id", "cell_type_mmc_raw"]],
    on="cell_id",
    how="left",
)

print("🔍 Annotations merged.")
print("Unique cell types:", adata.obs["cell_type_mmc_raw"].nunique())

print("🧪 Sanity check AFTER annotation merge:")
sanity_check(adata, color="cell_type_mmc_raw")


# -----------------------------------------
# 4) Compute UMAP on concept embedding
# -----------------------------------------
neighbors_umap(adata, use_rep="concept_mean_embedding")


# -----------------------------------------
# 5) Plot and save (same pattern as your zhuang code)
# -----------------------------------------
color = "cell_type_mmc_raw"
title = f"Zeng UMAP (concept embedding) colored by {color}"

cwd = os.getcwd()
out_png = os.path.join(cwd, f"zeng_umap_concept_{color}.png")
# If you don't want to save, set out_png = None
# out_png = None

plot(
    adata,
    color=color,
    out_png=out_png,
    title=title,
    palette=cell_type_palette,
)

print(f"🎉 Done! Plot saved to: {out_png}")
