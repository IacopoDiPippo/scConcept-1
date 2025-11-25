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


def neighbors_umap(adata, use_rep="concept_mean_embedding"):
    emb = adata.obsm[use_rep]      # shape: (n_cells, 128)

    # 1) Run PCA on the embedding directly
    print(f"Running PCA on {use_rep} …")
    pca = sc.tl.pca(emb, n_comps=50, svd_solver="arpack", copy=True)
    
    # Store PCA result in adata.obsm
    adata.obsm["X_pca_emb"] = pca

    # 2) Neighbors + UMAP on PCA of embedding
    print("PP Neighbors on PCA of embedding…")
    sc.pp.neighbors(adata, use_rep="X_pca_emb")

    print("Computing UMAP…")
    sc.tl.umap(adata, min_dist=0.3, random_state=0)



    """print("PP Neighbors on PCA…")
    sc.pp.neighbors(adata, use_rep="X_pca")

    print("Computing UMAP…")
    sc.tl.umap(adata, min_dist=0.3, random_state=0)"""



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


cell_type_palette = {'ABCs': '#023fa5',
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
 'VLMCs': '#f79cd4'}

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

neighbors_umap(adata, use_rep="concept_mean_embedding") # use concept embeddingor or "concept_mean_embedding"

# ---------------------------------------------------
# Plot
# ---------------------------------------------------

title = "Zhuang UMAP (concept embedding)"
out_png = os.path.join(os.getcwd(), "zhuang_umap_concept_mean_pca.png")

plot_umap(
    adata,
    color="cell_type",
    title=title,
    palette=cell_type_palette,  # add your palette if you want
    out_png=out_png
)

print("🎉 DONE!")
