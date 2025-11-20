import pandas as pd
import mygene

# -----------------------------------
# 1. LOAD EXCEL AND EXTRACT GENE NAMES
# -----------------------------------

"""

df = pd.read_csv("datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv")

# Column should be named "Gene"

# Get the list of cell names (the first column)
cell_names = df.iloc[:, 0].tolist()

# In the new dataset, genes are the column names (except the first one)
genes = df.columns[1:].tolist()

print(f"Loaded {len(genes)} gene symbols.")"""
"""gene_col = [c for c in df.columns if c.lower().startswith("name")][0]

genes = df[gene_col].dropna().unique().tolist()
"""
import scanpy as sc

adata = sc.read_h5ad("/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/processed/ISD-1.h5ad")

# Extract gene symbols from var index
genes = adata.var.index.tolist()

print(f"Loaded {len(genes)} gene symbols from ISD-1.h5ad.")

# -----------------------------------
# 2. MAP TO ENSEMBL IDs (GENE-LEVEL)
# -----------------------------------

mg = mygene.MyGeneInfo()

results = mg.querymany(
    genes,
    scopes="symbol",
    fields="ensembl.gene",
    species="mouse"
)

mapped_rows = []

for entry in results:
    symbol = entry["query"]
    ensembl_id = None
    
    if "ensembl" in entry:
        ens = entry["ensembl"]
        
        # If multiple mappings, take the first gene ID
        if isinstance(ens, list):
            ids = [e["gene"] for e in ens if "gene" in e]
            if ids:
                ensembl_id = ids[0]
        else:
            ensembl_id = ens.get("gene")
    
    mapped_rows.append({
        "Ensembl_ID": ensembl_id
    })

mapped_df = pd.DataFrame(mapped_rows)

# -----------------------------------
# 3. SAVE OUTPUT
# -----------------------------------

OUTPUT = "Vizgen1000.csv"
mapped_df.to_csv(OUTPUT, index=False)

print(f"Saved mapped Ensembl IDs to {OUTPUT}")
