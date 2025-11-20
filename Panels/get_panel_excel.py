import pandas as pd
import mygene

# -----------------------------------
# 1. LOAD EXCEL AND EXTRACT GENE NAMES
# -----------------------------------

FILENAME = "eng  Allen Institute MERFISH whole-mouse-brain gene panels.xlsx"

df = pd.read_excel(FILENAME)

# Column should be named "Gene"
gene_col = [c for c in df.columns if c.lower().startswith("gene")][0]

genes = df[gene_col].dropna().unique().tolist()

print(f"Loaded {len(genes)} gene symbols.")

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
        "Gene_Symbol": symbol,
        "Ensembl_ID": ensembl_id
    })

mapped_df = pd.DataFrame(mapped_rows)

# -----------------------------------
# 3. SAVE OUTPUT
# -----------------------------------

OUTPUT = "Allen_Zeng_EnsemblIDs.csv"
mapped_df.to_csv(OUTPUT, index=False)

print(f"Saved mapped Ensembl IDs to {OUTPUT}")
