import pandas as pd
from mygene import MyGeneInfo

mg = MyGeneInfo()

def symbol_to_geneid(symbol):
    result = mg.query(symbol, scopes="symbol", fields="ensembl.gene", species="mouse")
    if "hits" not in result or len(result["hits"]) == 0:
        return None
    ens = result["hits"][0].get("ensembl")
    if isinstance(ens, dict):
        return ens.get("gene")
    if isinstance(ens, list):
        return ens[0].get("gene")
    return None

# Read your file
df = pd.read_csv("Vizgen_Gene List_mouse_brain_1000.csv")

# The gene symbol is in the SECOND column:
symbols = df["name"]

# Convert symbols → Ensembl gene IDs
gene_ids = symbols.apply(symbol_to_geneid)

# Save output
pd.DataFrame({"Ensembl_ID": gene_ids}).to_csv("Vizgen1000_geneIDs.csv", index=False)

print("DONE: Vizgen1000_geneIDs.csv created!")
