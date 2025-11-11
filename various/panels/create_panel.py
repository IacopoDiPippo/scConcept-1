import anndata as ad
import pandas as pd
from pathlib import Path

# ---- SETTINGS ----
adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zeng.h5ad" # change this
output_dir = Path("/mouse")      # your panels directory
panel_name = "ZengGenePanel.csv"       # name of the panel CSV to create

# ---- LOAD DATA ----
adata = ad.read_h5ad(adata_path)

# ---- GET ENSEMBL IDS ----
# Check where your Ensembl gene IDs are stored
if "gene_ids" in adata.var.columns:
    ensembl_ids = adata.var["gene_ids"].astype(str)
    print("✅ Using adata.var['gene_ids'] for Ensembl IDs.")
else:
    ensembl_ids = adata.var_names.astype(str)
    print("⚠️ Using adata.var_names for Ensembl IDs (no 'gene_ids' column found).")

# ---- CREATE DATAFRAME ----
panel_df = pd.DataFrame({"Ensembl_ID": ensembl_ids})

# ---- SAVE CSV ----
output_dir.mkdir(parents=True, exist_ok=True)
panel_path = output_dir / panel_name
panel_df.to_csv(panel_path, index=False)

print(f"✅ Saved panel with {len(panel_df)} genes to {panel_path}")
