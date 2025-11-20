import anndata as ad
import pandas as pd
from pathlib import Path

# ---- SETTINGS ----
adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/abc_atlas.h5ad"
output_dir = Path("/p/home/jusers/dipippo1/jureca/projects/scConcept-1/Panels/done_panels")       # where to save the panel
panel_name = "ZengGenePanel.csv"  # name of the CSV

print(f"🔍 Reading {adata_path} in backed mode (-r)...")
adata = ad.read_h5ad(adata_path, backed='r')

# ---- GET ENSEMBL IDS ----
if "gene_ids" in adata.var.columns:
    print("✅ Found 'gene_ids' column in .var, using that for Ensembl IDs.")
    ensembl_ids = adata.var["gene_ids"].astype(str)
else:
    print("⚠️ No 'gene_ids' column found, using .var_names instead.")
    ensembl_ids = adata.var_names.astype(str)

# ---- LOGGING PROGRESS ----
n_genes = len(ensembl_ids)
print(f"📊 Found {n_genes} genes in the dataset.")

# ---- CREATE PANEL DATAFRAME ----
panel_df = pd.DataFrame({"Ensembl_ID": ensembl_ids})

# ---- SAVE ----
output_dir.mkdir(parents=True, exist_ok=True)
panel_path = output_dir / panel_name
panel_df.to_csv(panel_path, index=False)

print(f"\n✅ Done! Saved panel with {len(panel_df)} genes to {panel_path}")
