import anndata as ad
import pandas as pd
from pathlib import Path

# ---- SETTINGS ----
adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/abc_atlas.h5ad.h5ad"
output_path = Path("/p/home/jusers/dipippo1/jureca/projects/scConcept-1/various/pc_gene_token_mapping.pkl")

print(f"🔍 Reading {adata_path} in backed mode (-r)...")
adata = ad.read_h5ad(adata_path, backed='r')

# ---- GET ENSEMBL IDS ----
if "gene_ids" in adata.var.columns:
    print("✅ Found 'gene_ids' column in .var, using that for Ensembl IDs.")
    genes = adata.var["gene_ids"].astype(str)
else:
    print("⚠️ No 'gene_ids' column found, using .var_names instead.")
    genes = adata.var_names.astype(str)

n_genes = len(genes)
print(f"📊 Found {n_genes} genes in total.")

# ---- CREATE TOKEN MAPPING ----
# Reserve 0 and 1 for special tokens
special_tokens = {"<cls>": 0, "<pad>": 1}

# Create mapping for genes starting from index 2
gene_to_token = {gene: i for i, gene in enumerate(genes, start=2)}
gene_to_token.update(special_tokens)

# Convert to Pandas Series (required by scConcept)
gene_series = pd.Series(gene_to_token)

# ---- SAVE FILE ----
output_path.parent.mkdir(parents=True, exist_ok=True)
gene_series.to_pickle(output_path)

print(f"✅ Saved gene-to-token mapping with {len(gene_series)} entries to {output_path}")

# ---- OPTIONAL CHECK ----
print("\n🔎 Sample entries:")
print(gene_series.head(10))
