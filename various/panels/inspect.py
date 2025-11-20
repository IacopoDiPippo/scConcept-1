import anndata as ad

adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/processed/ISD-1.h5ad"

# Load the AnnData object
adata = ad.read_h5ad(adata_path)

# --- Basic structure ---
print("=== AnnData Structure ===")
print(adata)

# --- Names and shapes ---
print("\n=== Names ===")
print("Observations (obs) columns:", adata.obs.columns.tolist())
print("Variables (var) columns:", adata.var.columns.tolist())

print("\n=== Shapes ===")
print(f"Number of cells (obs): {adata.n_obs}")
print(f"Number of genes/features (vars): {adata.n_vars}")

# --- Examples ---
# One example obs row
print("\n=== One example from obs ===")
print(adata.obs.iloc[0])

# One example var row
print("\n=== One example from var ===")
print(adata.var.iloc[0])
