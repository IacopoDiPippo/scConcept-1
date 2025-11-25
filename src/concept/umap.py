import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
def sanity_check(adata, color):
    print("\n======================")
    print("🔍   UMAP SANITY CHECK")
    print("======================")

    print(f"Number of cells: {adata.n_obs:,}")
    print(f"Number of features: {adata.n_vars:,}")

    # Count missing labels
    missing = adata.obs[color].isna().sum()
    print(f"Missing labels for '{color}': {missing:,}")

    # Show unique categories
    print("\nUnique categories found:")
    print(adata.obs[color].unique())

    # Check representation existence
    if "X_umap" in adata.obsm:
        print("\nUMAP coordinates shape:", adata.obsm["X_umap"].shape)
    else:
        print("\n⚠️  WARNING: No UMAP coordinates found (X_umap missing).")

    print("======================\n")


# -----------------------------------------
# UMAP & Plot helpers
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
            adata, color=color, title=title,
            frameon=False, palette=palette,
            show=False, return_fig=True
        )
    else:
        ax = sc.pl.umap(
            adata, color=color, title=title,
            frameon=False, show=False, return_fig=True
        )

    plt.show()

    if out_png:
        print(f"Saving figure to: {out_png}")
        ax.savefig(out_png, dpi=300, bbox_inches="tight")

    plt.close(ax)


# -----------------------------------------
# Palette
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
# Load your saved AnnData
# -----------------------------------------
adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zeng_concept_embeddings.h5ad"
print(f"📂 Loading: {adata_path}")

adata = sc.read(adata_path)

# -----------------------------------------
# 0. Sanity check BEFORE merging annotations
# -----------------------------------------
print("🧪 Sanity check BEFORE annotation merge:")
print("obs columns:", adata.obs.columns.tolist())
print("shape:", adata.shape)
print("obsm keys:", adata.obsm.keys())

# You must pick a column that already exists
# Example: show distribution of donor or subclass
if "subclass" in adata.obs:
    sanity_check(adata, color="subclass")
else:
    print("⚠️ 'subclass' not found in obs — skipping pre-merge plot.")


# -----------------------------------------
# Attach annotations
# (Change to your actual annotation file path)
# -----------------------------------------
annotations_file = "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv"
ann = pd.read_csv(annotations_file)



adata.obs["cell_id"] = adata.obs.index

adata.obs = adata.obs.merge(
    ann[["cell_id", "cell_type_mmc_raw"]],
    on="cell_id",
    how="left"
)

print("🔍 Annotations merged.")
print("Unique cell types:", adata.obs["cell_type_mmc_raw"].nunique())


# -----------------------------------------
# 2. Sanity check AFTER merging
# -----------------------------------------
print("🧪 Sanity check AFTER annotation merge (cell_type_mmc_raw):")

sanity_check(adata, color="cell_type_mmc_raw")
# -----------------------------------------
# Compute UMAP on concept embedding
# -----------------------------------------
neighbors_umap(adata, use_rep="concept_mean_embedding")


# -----------------------------------------
# Plotting
# -----------------------------------------
color = "cell_type_mmc_raw"
title = f"Zeng UMAP (concept embedding) colored by {color}"

# Save plot to current working directory
cwd = os.getcwd()
out_png = os.path.join(cwd, f"zeng_umap_concept_{color}.png")

plot(
    adata,
    color=color,
    out_png=out_png,
    title=title,
    palette=cell_type_palette
)

print(f"\n🎉 Plot saved at:\n{out_png}\n")
