import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------
# Helpers
# ---------------------------------------------------

def sanity_check(adata, color=None):
    print("\n======================")
    print("🔍   ANNDATA CHECK")
    print("======================")
    print(f"Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print("obs columns:", list(adata.obs.columns))
    print("obsm keys:", list(adata.obsm.keys()))

    if color is not None and color in adata.obs:
        print(f"\nColumn '{color}': {adata.obs[color].nunique()} unique values")
        print(f"Missing: {adata.obs[color].isna().sum():,}")
    elif color is not None:
        print(f"\n⚠️ Column '{color}' NOT FOUND")

    if "X_umap" in adata.obsm:
        print("UMAP present:", adata.obsm["X_umap"].shape)
    else:
        print("No UMAP yet")
    print("======================\n")


def neighbors_umap(adata, use_rep=None):
    print("PP Neighbors…")
    sc.pp.neighbors(adata, use_rep=use_rep)
    print("Computing UMAP…")
    sc.tl.umap(adata, min_dist=0.3, random_state=0)


def plot_umap(adata, color, palette=None, title=None, out_png=None):
    fig = sc.pl.umap(
        adata,
        color=color,
        title=title,
        frameon=False,
        palette=palette,
        show=False,
        return_fig=True,
    )
    plt.show()

    if out_png:
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        print(f"Saved: {out_png}")


# ---------------------------------------------------
# Load Zhuang embeddings
# ---------------------------------------------------

adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zhuang-ABCA-1_concept_embeddings.h5ad"
print(f"📂 Loading {adata_path}")
adata = sc.read(adata_path)

sanity_check(adata)

# ---------------------------------------------------
# Prepare IDs
# ---------------------------------------------------

# The index of your adata is in .obs["cell_label"]
adata.obs["cell_id"] = adata.obs.index.astype(str)

# ---------------------------------------------------
# Subsample
# ---------------------------------------------------

print("📉 Subsampling to 100k…")
sc.pp.subsample(adata, n_obs=100000, random_state=0)
sanity_check(adata)

# ---------------------------------------------------
# Load and fix Zhuang annotations
# ---------------------------------------------------

annot_path = "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv"
print(f"📂 Loading annotation: {annot_path}")

ann = pd.read_csv(annot_path)

# Strip suffix "-Zhuang-ABCA1"
ann["base_id"] = ann["cell_id"].astype(str).str.split("-").str[0]

# Rename so merge works
ann = ann.rename(columns={"cell_type_mmc_raw": "cell_type"})

# ---------------------------------------------------
# Merge
# ---------------------------------------------------

adata.obs = adata.obs.merge(
    ann[["base_id", "cell_type"]],
    left_on="cell_id",
    right_on="base_id",
    how="left"
)

sanity_check(adata, color="cell_type")

# ---------------------------------------------------
# Compute UMAP on concept embedding
# ---------------------------------------------------

neighbors_umap(adata, use_rep="concept_mean_embedding")

# ---------------------------------------------------
# Plot
# ---------------------------------------------------

title = "Zhuang UMAP (concept embedding)"
out_png = os.path.join(os.getcwd(), "zhuang_umap_concept.png")

plot_umap(
    adata,
    color="cell_type",
    title=title,
    palette=None,  # add your palette if you want
    out_png=out_png
)

print("🎉 DONE!")
