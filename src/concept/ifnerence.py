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
from concept.data.datamodules import InferenceDataset, make_inference_dataloader
import os
import torch
import lightning as L
import pandas as pd
import anndata as ad

from hydra import compose, initialize
from omegaconf import OmegaConf

from lamin_dataloader.dataset import GeneIdTokenizer

# IMPORT CORRETTO DEL MODELLO - aggiorna il path a quello reale
from concept.models.bi_encoder_contrastive import BiEncoderContrastiveModel  # <-- METTI QUI IL TUO PATH

from inference_dataset import make_inference_dataloader  # oppure importa le classi se le hai nello stesso file


def main():
    # -------------------------
    # 1) Hydra config
    # -------------------------
    print("Loading Hydra config...")
    with initialize(version_base=None, config_path="conf"):
        cfg = compose(config_name="config")

    # -------------------------
    # 2) Gene mapping + tokenizer
    # -------------------------
    gene_mapping = pd.read_pickle(cfg.PATH.gene_mapping_path).to_dict()
    tokenizer = GeneIdTokenizer(gene_mapping)

    # -------------------------
    # 3) Path del file Zeng.h5ad
    # -------------------------
    val_files = cfg.PATH.SPLIT["val"]
    if isinstance(val_files, str):
        val_files = [val_files]
    assert len(val_files) == 1, f"Mi aspetto un solo file di validazione, trovato: {val_files}"

    h5ad_path = os.path.join(cfg.PATH.ADATA_PATH, val_files[0])
    print(f"Inference on: {h5ad_path}")

    # -------------------------
    # 4) Dataloader di inference
    # -------------------------
    dataloader, adata = make_inference_dataloader(
        h5ad_path=h5ad_path,
        tokenizer=tokenizer,
        batch_size=256,
        num_workers=4,
    )

    # -------------------------
    # 5) Carica il modello
    # -------------------------
    ckpt_path = os.path.join(
        cfg.PATH.CHECKPOINT_ROOT,
        "xtpijc3k",   # il run id che hai detto
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

    # -------------------------
    # 6) Predict
    # -------------------------
    trainer = L.Trainer(accelerator="gpu", devices=1, logger=False)
    print("\nRunning prediction...")
    preds = trainer.predict(model, dataloaders=dataloader)

    print("Collecting embeddings...")
    cls_emb = torch.cat([p["cls_cell_emb"] for p in preds], dim=0).cpu().numpy()
    mean_emb = torch.cat([p["mean_cell_emb"] for p in preds], dim=0).cpu().numpy()

    assert adata.n_obs == cls_emb.shape[0], (
        f"Num cells mismatch: adata.n_obs={adata.n_obs}, cls_emb={cls_emb.shape[0]}"
    )

    # -------------------------
    # 7) Salva su AnnData
    # -------------------------
    adata.obsm["concept_cls_embedding"] = cls_emb
    adata.obsm["concept_mean_embedding"] = mean_emb

    out_path = h5ad_path.replace(".h5ad", "_concept_embeddings.h5ad")
    adata.write(out_path)

    print(f"\n✅ Saved embeddings to:\n{out_path}\n")


if __name__ == "__main__":
    main()
