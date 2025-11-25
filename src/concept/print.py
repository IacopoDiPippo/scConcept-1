import scanpy as sc
adata = sc.read("/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/raw/Zhuang-ABCA-1_concept_embeddings.h5ad")

print("First 20 cell IDs:")
print(adata.obs.index[:20])
import pandas as pd

ann = pd.read_csv("/p/project1/hai_fzj_bda/salg1/cellseg-benchmark/data_dir/samples/zhuang/results/merfish/cell_type_annotation/adata_obs_annotated.csv")

print("Annotation cell IDs:")
print(ann["cell_id"][:20])
