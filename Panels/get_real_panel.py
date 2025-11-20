import pandas as pd
import ast

# ==========================
# 1. LOAD YOUR FILE
# ==========================

FILENAME = "gene.csv"   # <-- update this
df = pd.read_csv(FILENAME)

# ==========================
# 2. PARSE ID COLUMN SAFELY
# ==========================

all_ids = []

for val in df["gene_identifier"]:
    if pd.isna(val):
        continue
    
    # CASE A: It's a single ID (normal string)
    if isinstance(val, str) and not val.startswith("["):
        all_ids.append(val)
        continue
    
    # CASE B: It's a list stored as a string => convert using ast.literal_eval
    try:
        parsed = ast.literal_eval(val)
        if isinstance(parsed, list):
            for item in parsed:
                if isinstance(item, str):
                    all_ids.append(item)
    except:
        # fallback: add raw value
        all_ids.append(str(val))

# ==========================
# 3. REMOVE DUPLICATES
# ==========================

unique_ids = list(dict.fromkeys(all_ids))  # preserves order

# ==========================
# 4. SAVE OUTPUT
# ==========================

out_df = pd.DataFrame({"Ensembl_ID": unique_ids})
OUTPUT = "Ensembl_IDs.csv"
out_df.to_csv(OUTPUT, index=False)

print(f"Extracted {len(unique_ids)} unique Ensembl IDs → saved to {OUTPUT}")
