import anndata as ad
import pandas as pd
from joblib import dump
from pathlib import Path
from sklearn.preprocessing import LabelEncoder

# ---- SETTINGS ----
adata_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/abc_atlas.h5ad"
output_path = Path("/p/home/jusers/dipippo1/jureca/projects/scConcept-1/various/label_encoders.joblib")

print(f"🔍 Reading {adata_path}...")
adata = ad.read_h5ad(adata_path, backed='r')

# ---- DEFINE COLUMNS ----
cols_to_encode = []
if "cell_type" in adata.obs.columns:
    cols_to_encode.append("cell_type")
else:
    print("⚠️ No 'cell_type' found in adata.obs")

# some datasets use 'donor_label' or 'donor_id'
for donor_col in ["donor_id", "donor_label"]:
    if donor_col in adata.obs.columns:
        cols_to_encode.append(donor_col)
        break
else:
    print("⚠️ No donor ID column found ('donor_id' or 'donor_label')")

if not cols_to_encode:
    raise ValueError("❌ No valid categorical columns found for label encoding!")

print(f"📊 Will encode: {cols_to_encode}")

# ---- BUILD LABEL ENCODERS ----
encoders = {}
for col in cols_to_encode:
    le = LabelEncoder()
    values = adata.obs[col].astype(str).fillna("Unknown")
    le.fit(values)
    encoders[col] = dict(zip(le.classes_, le.transform(le.classes_)))
    print(f"✅ Encoded '{col}' with {len(le.classes_)} unique classes.")

# ---- SAVE ----
output_path.parent.mkdir(parents=True, exist_ok=True)
dump(encoders, output_path)
print(f"💾 Saved label encoders to {output_path}")
