import anndata as ad
import pandas as pd

adata = ad.read_h5ad("/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/abc_atlas.h5ad")
gene_mapping = pd.read_pickle(/p/scratch/cjinm16/dipippo1/scConcept/token_mapping/pc_gene_token_mapping.pkl).to_dict()

common = set(adata.var_names).intersection(set(gene_mapping.keys()))
print("adata genes:", adata.n_vars)
print("mapping genes:", len(gene_mapping))
print("common:", len(common))