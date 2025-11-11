import os
from pathlib import Path
import scanpy as sc
from concept.scConcept import scConcept

# ------------------------------------------------------------------------------
# 1. Setup cache directory
# ------------------------------------------------------------------------------
cache_dir = Path("./cache/")
os.makedirs(cache_dir, exist_ok=True)

# ------------------------------------------------------------------------------
# 2. Download dataset if not already present
# ------------------------------------------------------------------------------
filename = cache_dir / "cite_gex_processed_training.h5ad"
url = "https://openproblems-bio.s3.amazonaws.com/public/explore/cite/cite_gex_processed_training.h5ad"

if not os.path.exists(filename):
    import urllib.request
    print(f"Downloading {filename} ...")
    urllib.request.urlretrieve(url, filename)
else:
    print(f"{filename} already exists, skipping download.")

# ------------------------------------------------------------------------------
# 3. Load dataset
# ------------------------------------------------------------------------------
adata = sc.read(filename)
print("Loaded dataset:")
print(adata)

# ------------------------------------------------------------------------------
# 4. Load pretrained scConcept model
# ------------------------------------------------------------------------------
concept = scConcept(cache_dir=cache_dir)
concept.load_config_and_model(model_name="Corpus-30M")

# ------------------------------------------------------------------------------
# 5. Extract embeddings
# ------------------------------------------------------------------------------
result = concept.extract_embeddings(
    adata=adata,
    batch_size=32,
    gene_id_column="gene_ids",
)

print(f"CLS embeddings: {result['cls_cell_emb'].shape}")
print(f"Mean embeddings: {result['mean_cell_emb'].shape}")

# Add embeddings to AnnData object
adata.obsm["X_scConcept"] = result["cls_cell_emb"]

# ------------------------------------------------------------------------------
# 6. Compute UMAP (optional visualization)
# ------------------------------------------------------------------------------
print("Computing neighbors and UMAP...")
sc.pp.neighbors(adata, use_rep="X_scConcept")
sc.tl.umap(adata)
sc.pl.umap(adata, color="cell_type", show=True)

# ------------------------------------------------------------------------------
# 7. Save embeddings for later use
# ------------------------------------------------------------------------------
output_file = cache_dir / "adata_with_scConcept_embeddings.h5ad"
adata.write(output_file)
print(f"Embeddings saved to: {output_file}")
