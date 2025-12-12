import os
import os
import torch
import lightning as L
import anndata as ad
import pandas as pd
from hydra import compose, initialize
from omegaconf import OmegaConf

from model import BiEncoderContrastiveModel
from lamin_dataloader.dataset import GeneIdTokenizer
import os
import torch
import lightning as L
import pandas as pd
import anndata as ad

from hydra import compose, initialize
from omegaconf import OmegaConf

from lamin_dataloader.dataset import GeneIdTokenizer

# IMPORT CORRETTO DEL MODELLO - aggiorna il path a quello reale
from model import BiEncoderContrastiveModel  # <-- METTI QUI IL TUO PATH



from lamin_dataloader.dataset import GeneIdTokenizer

from data.datamodules import AnnDataModule, InferenceAnnDataModule


def main():
    torch.set_float32_matmul_precision("high")
    # ----------------------------
    # 1. Load Hydra config
    # ----------------------------
    print("Loading config...")
    with initialize(version_base=None, config_path="conf"):
        cfg = compose(config_name="config")

    print("Original config loaded.")
    print(OmegaConf.to_yaml(cfg))


    # ----------------------------
    # 3. Load tokenizer
    # ----------------------------
    gene_mapping = pd.read_pickle(cfg.PATH.gene_mapping_path).to_dict()
    tokenizer = GeneIdTokenizer(gene_mapping)

    # ----------------------------
    # 4. Path to Zeng.h5ad
    # ----------------------------
    val_files = cfg.PATH.SPLIT["val"]
    if isinstance(val_files, str):
        val_files = [val_files]
    assert len(val_files) == 1

    h5ad_path = os.path.join(cfg.PATH.ADATA_PATH, val_files[0])
    h5ad_path = "/p/project1/hai_fzj_bda/spitzer2/point_transformer/data/processed/ISD-1.h5ad"
    print(f"Inference on: {h5ad_path}")
    cfg.datamodule.columns = []
    # ----------------------------
    # 4b. PREPROCESS: convert gene symbols → Ensembl IDs
    # ----------------------------
    print("Loading AnnData for preprocessing...")
    adata_raw = ad.read_h5ad(h5ad_path)
    
    print(f"Found {adata_raw.n_vars} gene symbols.")
    
    # ----------------------------
    # Load LOCAL precomputed mapping
    # ----------------------------
    mapping_path = "/p/home/jusers/dipippo1/jureca/projects/scConcept-1/Panels/done_panels/ISD2.csv"
    print(f"Loading local Ensembl mapping: {mapping_path}")
    
    ens_map = pd.read_csv(mapping_path)
    
    if len(ens_map) != adata_raw.n_vars:
        raise ValueError("Mapping rows do NOT match number of genes in AnnData!")
    
    # Assign Ensembl IDs
    adata_raw.var["ensembl_id"] = ens_map["Ensembl_ID"].values
    adata_raw.var_names = adata_raw.var["ensembl_id"].astype(str)  # Set index
    
    # Remove unmapped
    missing = (adata_raw.var_names == "nan").sum() + adata_raw.var_names.isna().sum()
    if missing > 0:
        print(f"⚠️ {missing} genes unmapped — removing them.")
        keep = ~(adata_raw.var_names.isna() | (adata_raw.var_names == "nan"))
        adata_raw = adata_raw[:, keep].copy()
    
    print(f"Final genes after mapping: {adata_raw.n_vars}")
    
    # Remove duplicate column to satisfy H5AD writing rules
    if "ensembl_id" in adata_raw.var.columns:
        del adata_raw.var["ensembl_id"]
    
    mapped_h5ad_path = h5ad_path.replace(".h5ad", "_ENSEMBL.h5ad")
    adata_raw.write(mapped_h5ad_path)
    
    print(f"Saved mapped AnnData to: {mapped_h5ad_path}")

    


    # ----------------------------
    # 5. Build inference DataModule
    # ----------------------------
    print("Building inference datamodule...")

    infer_dm = InferenceAnnDataModule(
        adata_path=mapped_h5ad_path,
        tokenizer=tokenizer,
        columns=cfg.datamodule.columns,
        normalization=cfg.datamodule.normalization,
        max_tokens=cfg.datamodule.dataset.train.max_tokens,
        batch_size=256,
        num_workers=4,
    )

    infer_dm.setup()
    val_loader = infer_dm.predict_dataloader()

    # ----------------------------
    # 6. Load trained model
    # ----------------------------
    ckpt_path = os.path.join(
        cfg.PATH.CHECKPOINT_ROOT,
        "xtpijc3k",
        "epochs",
        "last.ckpt",
    )
    print(f"Loading checkpoint: {ckpt_path}")

    model = BiEncoderContrastiveModel.load_from_checkpoint(
        ckpt_path,
        config=cfg.model,
        pad_token_id=gene_mapping["<pad>"],
        cls_token_id=gene_mapping["<cls>"],
        vocab_size=len(gene_mapping),
        precomp_embs_key=None,
        world_size=1,
        val_loader_names=[],
        strict=False,
    )
    model.eval()

    # ----------------------------
    # 7. Predict embeddings
    # ----------------------------
    trainer = L.Trainer(accelerator="gpu", devices=1, precision="16-mixed", logger=False)
    print("\nRunning prediction...")

    preds = trainer.predict(model, dataloaders=val_loader)

    cls_emb = torch.cat([p["cls_cell_emb"] for p in preds], dim=0).cpu().numpy()
    mean_emb = torch.cat([p["mean_cell_emb"] for p in preds], dim=0).cpu().numpy()

    # ----------------------------
    # 8. Save to AnnData
    # ----------------------------
    adata = ad.read_h5ad(h5ad_path)

    adata.obsm["concept_cls_embedding"] = cls_emb
    adata.obsm["concept_mean_embedding"] = mean_emb

    out_path = h5ad_path.replace(".h5ad", "_concept_embeddings.h5ad")
    adata.write(out_path)

    print(f"\n🎉 Saved to: {out_path}\n")


if __name__ == "__main__":
    main()
