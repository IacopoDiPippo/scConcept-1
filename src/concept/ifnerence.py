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
from concept.data.datamodules import InferenceDataModule


# ----------------------------------------------------------
#                    MAIN INFERENCE SCRIPT
# ----------------------------------------------------------
def main():

    # 1) LOAD HYDRA CONFIG
    print("Loading Hydra config...")
    with initialize(version_base=None, config_path="conf"):
        cfg = compose(config_name="config")

    # 2) LOAD TOKENIZER
    gene_mapping = pd.read_pickle(cfg.PATH.gene_mapping_path).to_dict()
    tokenizer = GeneIdTokenizer(gene_mapping)

    # choose validation file (Zeng.h5ad)
    val_file = cfg.PATH.SPLIT["val"]
    if isinstance(val_file, list):
        assert len(val_file) == 1, "We expect only one validation file for inference."
        val_file = val_file[0]

    h5ad_path = os.path.join(cfg.PATH.ADATA_PATH, val_file)
    print(f"Inference on: {h5ad_path}")

    # 3) BUILD INFERENCE DATAMODULE
    dm = InferenceDataModule(
        h5ad_path=h5ad_path,
        tokenizer=tokenizer,
        batch_size=256
    )
    dm.setup()

    # 4) LOAD MODEL CHECKPOINT
    ckpt_path = os.path.join(
        cfg.PATH.CHECKPOINT_ROOT,
        "xtpijc3k",
        "epochs",
        "last.ckpt"
    )
    print(f"Loading model checkpoint:\n{ckpt_path}")

    model = BiEncoderContrastiveModel.load_from_checkpoint(
        ckpt_path,
        config=cfg.model,
        pad_token_id=gene_mapping["<pad>"],
        cls_token_id=gene_mapping["<cls>"],
        vocab_size=len(gene_mapping),
        precomp_embs_key=None,
        world_size=1,
        val_loader_names=[],
        strict=False
    )
    model.eval()

    # 5) RUN PREDICT
    trainer = L.Trainer(accelerator="gpu", devices=1, logger=False)
    preds = trainer.predict(model, datamodule=dm)

    # COLLECT EMBEDDINGS
    cls_emb = torch.cat([p["cls_cell_emb"] for p in preds], dim=0).cpu().numpy()
    mean_emb = torch.cat([p["mean_cell_emb"] for p in preds], dim=0).cpu().numpy()

    # 6) SAVE OUTPUT
    adata = dm.adata
    adata.obsm["concept_cls_embedding"] = cls_emb
    adata.obsm["concept_mean_embedding"] = mean_emb

    out_path = h5ad_path.replace(".h5ad", "_concept_embeddings.h5ad")
    adata.write(out_path)

    print(f"\n🎉 Saved embeddings to:\n{out_path}\n")


if __name__ == "__main__":
    main()