import argparse
import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import accuracy_score, f1_score
from sklearn.preprocessing import LabelEncoder

from torch.utils.data import DataLoader, TensorDataset
import pytorch_lightning as pl
import torch 
import torch.nn as nn 
from torchvision.ops import MLP
from torchmetrics.classification import MulticlassAccuracy, MulticlassF1Score
from sklearn.model_selection import train_test_split

class MLPClassifier(pl.LightningModule):
    def __init__(self, 
                 input_size, 
                 n_classes , 
                 hidden_units=[128, 128], 
                 class_weights=None, 
                 lr = 1e-3,
                 weight_decay=0.01,
                 dropout=0.0):
        super().__init__()
        
        self.lr = lr
        self.weight_decay = weight_decay
        
        self.train_acc = MulticlassF1Score(num_classes=n_classes, average='macro')
        self.val_acc = MulticlassF1Score(num_classes=n_classes, average='macro')
        self.test_acc = MulticlassF1Score(num_classes=n_classes, average='macro')
        
        if class_weights is None:
            self.class_weights = torch.ones(n_classes, device=self.device, dtype=torch.float32)
        else:
            self.class_weights = torch.tensor(class_weights, device=self.device, dtype=torch.float32)
        
        if hidden_units is None or len(hidden_units) == 0:
            self.model = nn.Sequential(
                nn.Dropout(dropout), 
                nn.Linear(input_size, n_classes)
                )
        else:
            self.model = nn.Sequential(
                MLP(input_size, hidden_units, dropout=dropout), 
                nn.Linear(hidden_units[-1], n_classes)
                )

    def forward(self, x):
        x = self.model(x)
        return x

    def training_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x)
        loss = nn.functional.cross_entropy(logits, y, weight=self.class_weights.to(self.device))
        preds = torch.argmax(logits, dim=1)
        self.train_acc(preds, y)
        self.log("train_loss", loss, prog_bar=True)
        self.log("train_f1_macro", self.train_acc, prog_bar=True, on_step=True, on_epoch=True)
        return loss

    
    def validation_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x)
        loss = nn.functional.cross_entropy(logits, y, weight=self.class_weights.to(self.device))
        preds = torch.argmax(logits, dim=1)
        self.val_acc(preds, y)
        self.log("val_loss", loss, prog_bar=True)
        self.log("val_f1_macro", self.val_acc, prog_bar=True, on_step=False, on_epoch=True)
        return loss

    def test_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x)
        loss = nn.functional.cross_entropy(logits, y, weight=self.class_weights.to(self.device))
        preds = torch.argmax(logits, dim=1)
        self.test_acc(preds, y)
        self.log("test_loss", loss, prog_bar=True)
        self.log("test_f1_macro", self.test_acc, prog_bar=True, on_step=False, on_epoch=True)
        return loss
    
    # def predict_step(self, batch, batch_idx):
    #     x, y = batch
    #     logits = self(x)
    #     preds = torch.argmax(logits, dim=1)
    #     return preds

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        return optimizer



def _get_dataloader(embs, labels, batch_size=128, shuffle=False):
    dataset = TensorDataset(torch.tensor(embs).float(), torch.tensor(labels))
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)



class Classifier:
    def __init__(self, input_size, hidden_units=[], classes=None, use_class_weights=False):
        self.input_size = input_size
        self.hidden_units = hidden_units
        self.classes = list(classes)
        self.n_classes = len(classes)
        self.le = LabelEncoder()
        self.le.fit(classes)
        self.use_class_weights = use_class_weights


    def _get_model(self, **kwargs):
        return MLPClassifier(self.input_size, self.n_classes, hidden_units=self.hidden_units, **kwargs)
        
    def fit(self, train_embs, train_labels, val_embs, val_labels, max_epochs=None, min_epochs=0, batch_size=128, patience=5, min_delta=0.005, val_check_interval=500, **kwargs):

        if self.use_class_weights:
            class_weights = 1 / np.bincount(self.le.transform(train_labels))
            class_weights = class_weights / class_weights.max()
        else:
            class_weights = None
            
        train_dataloader = _get_dataloader(train_embs, self.le.transform(train_labels), batch_size=batch_size, shuffle=True)
        val_dataloader = _get_dataloader(val_embs, self.le.transform(val_labels), batch_size=batch_size, shuffle=False)

        self.model = self._get_model(class_weights=class_weights, **kwargs)

        early_stopping = pl.callbacks.EarlyStopping('val_f1_macro', mode='max', patience=patience, min_delta=min_delta, verbose=True)
        if torch.cuda.is_available():
            trainer = pl.Trainer(max_epochs=max_epochs, min_epochs=min_epochs, accelerator="gpu", devices=1, val_check_interval=val_check_interval, check_val_every_n_epoch=None, callbacks=[early_stopping], logger=False)
        else:
            trainer = pl.Trainer(max_epochs=max_epochs, min_epochs=min_epochs, val_check_interval=val_check_interval, check_val_every_n_epoch=None, callbacks=[early_stopping], logger=False)
            
        trainer.fit(self.model, train_dataloaders=train_dataloader, val_dataloaders=val_dataloader)
        self.epoch = trainer.current_epoch

    def predict(self, test_embs, pred_labels=None):
        self.model = self.model.eval()
        batch_size = 128
        test_dataloader = _get_dataloader(test_embs, np.zeros(test_embs.shape[0]), batch_size=batch_size, shuffle=False)
        if pred_labels is not None:
            pred_labels = self.le.transform(pred_labels)
        
        y_pred = []
        for x, _ in test_dataloader:
            x = x.to(self.model.device)
            logits = self.model(x)
            if pred_labels is not None:
                diff_val = logits.max() - logits.min() + 1
                logits[:, pred_labels] += diff_val
            y_pred.append(torch.argmax(logits, dim=1))
        y_pred = torch.cat(y_pred).cpu().numpy()
        return y_pred
    
    
    def save_model(self, path_dir):
        torch.save(self.model.state_dict(), os.path.join(path_dir, 'model.pth'))
        torch.save(self.le, os.path.join(path_dir, 'label_encoder.pth'))
        params = {
            'input_size': self.input_size,
            'hidden_units': self.hidden_units,
            'classes': self.classes,
            'use_class_weights': self.use_class_weights,
            'epoch': self.epoch
        }
        with open(os.path.join(path_dir, 'model_params.json'), 'w') as f:
            json.dump(params, f)
    
    
    @classmethod
    def from_checkpoint(cls, path_dir):
        with open(os.path.join(path_dir, 'model_params.json'), 'r') as f:
            params = json.load(f)
        
        epoch = params.pop('epoch')
        obj = cls(**params)
        obj.epoch = epoch
        obj.model = obj._get_model()
        obj.model.load_state_dict(torch.load(os.path.join(path_dir, 'model.pth')))
        obj.le = torch.load(os.path.join(path_dir, 'label_encoder.pth'))
        return obj