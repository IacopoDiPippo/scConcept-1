import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------
# Fixed palette exactly as in the notebook
# ------------------------------------------------------
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

# ------------------------------------------------------
# Helper functions
# ------------------------------------------------------
def sanity_check(adata, color=None):
    print("\n======================")
    print("🔍   ANNDATA CHECK")
    print("======================")
    print(f"Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print("obs columns:", list(adata.obs.columns))
    print("obsm keys:", list(adata.obsm.keys()))

    if color and color in adata.obs:
        miss = adata.obs[color].isna().sum()
        print(f"\nColumn '{color}': {adata.obs[color].nunique()} unique, {miss:,} missing")
    elif color:
        print(f"⚠️ Column '{color}' not found!")

    if "X_umap" in adata.obsm:
        print("UMAP already exists →", adata.obsm["X_umap"].shape)
    else:
        print("No UMAP yet")
    print("======================\n")


def neighbors_umap(adata, use_rep=None):
    print("PP Neighbors")
    sc.pp.neighbors(adata, use_rep=use_rep)

    print("UMAP")
    sc.tl.umap(adata, min_dist=0.3, random_state=0)


def plot_umap(adata, color, palette, title, out_png):
    fig = sc.pl.umap(
        adata,
        color=color,
        palette=palette,
        title=title,
        frameon=False,
        show=False,
        return_fig=True,
    )
    fig.show()

    if out_png:
        print(f"Saving figure to: {out_png}")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------
# MAIN SCRIPT
# ------------------------------------------------------
def main():

    # 1) Load Zhuang WITH concept embeddings
    adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zhuang-ABCA-1_concept_embeddings.h5ad"
    print(f"📂 Loading: {adata_path}")
    adata = sc.read(adata_path)

    print("🧪 Sanity check BEFORE:")
    sanity_check(adata)

    # 2) Add cell_id (needed for annotation merge)
    adata.obs["cell_id"] = adata.obs.index

    # 3) Subsample to 100k (like notebook)
    print("📉 Subsampling to 100,000 cells")
    sc.pp.subsample(adata, n_obs=100000, random_state=0)

    print("🧪 AFTER subsample:")
    sanity_check(adata)

    # 4) Merge Zhuang MERFISH annotations
    annot_file = "/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv"
    print(f"📂 Loading annotation: {annot_file}")
    ann = pd.read_csv(annot_file)

    adata.obs = adata.obs.merge(
        ann[["cell_id", "cell_type_mmc_raw"]],
        on="cell_id",
        how="left",
    )

    print("🔍 Annotations merged.")
    print("Unique cell types:", adata.obs["cell_type_mmc_raw"].nunique())

    print("🧪 AFTER annotations:")
    sanity_check(adata, color="cell_type_mmc_raw")

    # 5) Compute neighbors/UMAP on *concept_mean_embedding*
    neighbors_umap(adata, use_rep="concept_mean_embedding")

    print("🧪 AFTER UMAP:")
    sanity_check(adata, color="cell_type_mmc_raw")

    # 6) Plot UMAP
    title = "Zhuang UMAP (concept_mean_embedding) — cell_type_mmc_raw"
    out_png = os.path.join(os.getcwd(), "zhuang_umap_concept.png")

    plot_umap(
        adata,
        color="cell_type_mmc_raw",
        palette=cell_type_palette,
        title=title,
        out_png=out_png,
    )

    print(f"\n🎉 DONE — figure saved to: {out_png}\n")


if __name__ == "__main__":
    main()
