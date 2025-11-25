import os
import torch
import lightning as L
import anndata as ad
import pandas as pd


from model import BiEncoderContrastiveModel
from lamin_dataloader.dataset import GeneIdTokenizer
from concept.data.datamodules import InferenceDataModule


def main():

    # ----------------------------------------------
    # Load gene mapping + tokenizer
    # ----------------------------------------------
    gene_mapping = pd.read_pickle("/path/to/gene_mapping.pkl").to_dict()
    tokenizer = GeneIdTokenizer(gene_mapping)

    h5ad_path = "/p/project1/.../Zeng.h5ad"

    # ----------------------------------------------
    # Build inference datamodule
    # ----------------------------------------------
    dm = InferenceDataModule(h5ad_path, tokenizer, batch_size=256)
    dm.setup()

    # ----------------------------------------------
    # Load model checkpoint
    # ----------------------------------------------
    ckpt_path = "/p/scratch/.../xtpijc3k/epochs/last.ckpt"

    model = BiEncoderContrastiveModel.load_from_checkpoint(
        ckpt_path,
        config=cfg.model,
        pad_token_id=gene_mapping["<pad>"],
        cls_token_id=gene_mapping["<cls>"],
        vocab_size=len(gene_mapping),
        strict=False
    )
    model.eval()

    # ----------------------------------------------
    # Predict
    # ----------------------------------------------
    trainer = L.Trainer(accelerator="gpu", devices=1, logger=False)
    preds = trainer.predict(model, datamodule=dm)

    cls_emb = torch.cat([p["cls_cell_emb"] for p in preds]).cpu().numpy()
    mean_emb = torch.cat([p["mean_cell_emb"] for p in preds]).cpu().numpy()

    # ----------------------------------------------
    # Save new AnnData
    # ----------------------------------------------
    adata = dm.adata
    adata.obsm["concept_cls_embedding"] = cls_emb
    adata.obsm["concept_mean_embedding"] = mean_emb

    out_path = h5ad_path.replace(".h5ad", "_concept_embeddings.h5ad")
    adata.write(out_path)

    print(f"Saved embeddings to:\n{out_path}")


if __name__ == "__main__":
    main()
