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
    print(f"Inference on: {h5ad_path}")

    # ----------------------------
    # 5. Build inference DataModule
    # ----------------------------
    print("Building inference datamodule...")

    infer_dm = InferenceAnnDataModule(
        adata_path=h5ad_path,
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
    trainer = L.Trainer(accelerator="gpu", devices=1, logger=False)
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
