import argparse
import numpy as np
import scanpy as sc
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import train_test_split

from ct_rep.utils.sklearn_classifiers import knn_on_embeddings
from ct_rep.utils.pytorch_mlp import Classifier
import os
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import torch
import torch.nn.functional as F
from ct_rep.utils.utils import split_train_test

def evaluate(train_embs, train_labels, test_embs, test_labels, classifier, metric, val_frac=0.1, pred_within_test_labels=False, lr=1e-3, n_neighbors=10):
    assert classifier in ['knn', 'linear', 'mlp'], 'classifier must be one of knn, linear, mlp'
    
    print(f'train split: {len(train_embs)} cells containing {train_labels.nunique()} cell types, embeddings of size {train_embs.shape[1]}')
    print(f'test split: {len(test_embs)} cells containing {test_labels.nunique()} cell types, embeddings of size {test_embs.shape[1]}')
    
    train_embs, val_embs, train_labels, val_labels = train_test_split(train_embs, train_labels, test_size=val_frac, stratify=train_labels, random_state=42)

    # KNN
    if classifier == 'knn':
        y_pred = knn_on_embeddings(
                            train_embs,
                            test_embs,
                            train_labels,
                            test_labels,
                            n_neighbors = n_neighbors, 
                            metric = metric,
                            pred_labels=np.unique(test_labels) if pred_within_test_labels else None
                            )
    else:
        if classifier == 'linear':
            hidden_units = []
            dropout=0.1
            weight_decay=0.01
        elif classifier == 'mlp':
            hidden_units = [512, 256]
            dropout=0.2
            weight_decay=0.02
            
        if train_embs.shape[1] > 1024 or len(train_labels) < 5000:
            print('Training with high regularization')
            dropout = dropout * 2
            weight_decay = weight_decay * 2

        clf = Classifier(input_size=train_embs.shape[1],
                        hidden_units=hidden_units,
                        classes=np.unique(train_labels),
                        use_class_weights=False,
                        )
        clf.fit(train_embs, 
                train_labels, 
                val_embs=val_embs,
                val_labels=val_labels,
                max_epochs=1000,
                min_epochs=0,
                batch_size=256,
                patience=5,
                min_delta=0.003,
                val_check_interval=300,
                lr=lr,
                dropout=dropout, 
                weight_decay=weight_decay
                )
        y_pred = clf.predict(test_embs,
                             pred_labels=np.unique(test_labels) if pred_within_test_labels else None,
                             )
        y_pred = clf.le.inverse_transform(y_pred)

    accuracy = accuracy_score(test_labels, y_pred)
    f1_score_weighted = f1_score(test_labels, y_pred, labels=np.unique(test_labels), average='weighted')
    f1_score_macro = f1_score(test_labels, y_pred, labels=np.unique(test_labels), average='macro')
    print(f'{classifier} Accuracy: {accuracy:.3f}')
    print(f'{classifier} F1-Weighted: {f1_score_weighted:.3f}')
    print(f'{classifier} F1-Macro: {f1_score_macro:.3f}')
    
    return accuracy, f1_score_weighted, f1_score_macro, y_pred


def run(adata_train,
        adata_test,
        split_key_ref, 
        split_key_query, 
        split_value_ref, 
        split_value_query, 
        label_key_ref, 
        label_key_query, 
        classifier,
        metric='cosine',
        train_size=None, 
        lr=1e-3, 
        min_count=50, 
        write_report=False
        ):

    
    # remove low frequency cell types
    freq = adata_train.obs[label_key_ref].value_counts(normalize=False)
    print(f'removing {list(freq[freq < min_count].index)} from dataset')
    adata_train = adata_train[adata_train.obs[label_key_ref].isin(freq[freq > min_count].index)]
    adata_test = adata_test[adata_test.obs[label_key_query].isin(freq[freq > min_count].index)]
        
    adata_train, adata_test = split_train_test(adata_train, adata_test, split_key_ref, split_key_query, split_value_ref, split_value_query)
    
    adata_train = adata_train.copy()
    adata_test = adata_test.copy()
    if metric == 'cosine':
        adata_train.obsm['X_emb'] = F.normalize(torch.tensor(adata_train.obsm['X_emb']), p=2, dim=1).numpy()
        adata_test.obsm['X_emb'] = F.normalize(torch.tensor(adata_test.obsm['X_emb']), p=2, dim=1).numpy()
        
    train_embs = adata_train.obsm['X_emb']
    test_embs = adata_test.obsm['X_emb']
    train_labels = adata_train.obs[label_key_ref]
    test_labels = adata_test.obs[label_key_query]
    
    if train_size is not None:
        # adata_balanced = balance_anndata(adata_train, min_cells_per_class, label_key)
        sc.pp.subsample(adata_train, n_obs=train_size, copy=False, random_state=42)
        # adata_train = ad.concat([adata_train, adata_balanced], axis=0)

    print(train_embs)
    accuracy, f1_score_weighted, f1_score_macro, y_pred = evaluate(train_embs, train_labels, test_embs, test_labels, classifier, metric=metric, pred_within_test_labels=False, lr=lr)
    print(f'{accuracy:.3f}, {f1_score_weighted:.3f}, {f1_score_macro:.3f}')
    
    # report = classification_report(test_labels, y_pred, labels=np.unique(test_labels))
    
    # if write_report and adata_path_ref == adata_path_query:
    #     # print(report)
    #     report_path = os.path.join(os.path.dirname(emb_path_ref), f'classification_report_{classifier}_{split_key_ref}.txt')
    #     if not os.path.exists(report_path):
    #         with open(report_path, "w") as f:
    #             f.write(report)
    #         print(f"Report written to {report_path}")
    #     else:
    #         raise ValueError(f"Path {report_path} already exists!")

    #     classes_ = np.unique(test_labels)
    #     cm = confusion_matrix(test_labels, y_pred, labels=classes_)
    #     fig, ax = plt.subplots(figsize=(len(classes_)//1.5, len(classes_)//1.5))
    #     disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=classes_)
    #     disp.plot(xticks_rotation='vertical', ax=ax)
    #     plt.savefig(os.path.join(os.path.dirname(emb_path_ref), f'confusion_matrix_{classifier}_{split_key_ref}.png'), bbox_inches="tight")
        
    # return accuracy, f1_score_weighted, f1_score_macro

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--classifier', type=str, default='knn', help='classifier to use for classification')
    parser.add_argument("--adata_path_ref", type=str, required=True, help="adata path ref")
    parser.add_argument("--adata_path_query", type=str, required=True, help="adata path query")
    parser.add_argument("--emb_path_ref", type=str, required=True, help="embedding path")
    parser.add_argument("--emb_path_query", type=str, required=True, help="embedding path")
    parser.add_argument('--split_key_ref', type=str, default='batch', help='The key from adata.obs to use for splitting')
    parser.add_argument('--split_key_query', type=str, default='batch', help='The key from adata.obs to use for splitting')
    parser.add_argument('--split_value_ref', type=str, default='train', help='split to train over')
    parser.add_argument('--split_value_query', type=str, default='test', help='split to test over')
    parser.add_argument("--label_key_ref", type=str, default='cell_type', help="key to use for labels")
    parser.add_argument("--label_key_query", type=str, default='cell_type', help="key to use for labels")
    parser.add_argument("--metric", type=str, default='cosine', help="metric to use for knn")
    parser.add_argument('--train_size', type=str, default=None, help='number of training samples to use for training')
    parser.add_argument("--lr", type=float, default=1e-3, help="learning rate")
    parser.add_argument("--min_count", type=int, default=50, help="minimum number of cells per cell type")
    parser.add_argument("--write_report", action='store_true', default=False, help="write classification report")
    args = parser.parse_args()
    
    if os.path.exists(args.adata_path_ref) and os.path.exists(args.emb_path_ref) and os.path.exists(args.adata_path_query) and os.path.exists(args.emb_path_query):
        adata_train = sc.read_h5ad(args.adata_path_ref)
        adata_test = sc.read_h5ad(args.adata_path_query)
        adata_train.obsm['X_emb'] = np.load(args.emb_path_ref)
        adata_test.obsm['X_emb'] = np.load(args.emb_path_query)
        accuracy, f1_score_weighted, f1_score_macro = run(
            adata_train,
            adata_test,
            args.split_key_ref,
            args.split_key_query,
            args.split_value_ref,
            args.split_value_query,
            args.label_key_ref,
            args.label_key_query,
            metric=args.metric,
            train_size=args.train_size,
            classifier=args.classifier,
            lr=args.lr,
            min_count=args.min_count,
            # write_report=args.write_report
        )
    else:
        print(f'One of adata_path_ref, adata_path_query, emb_path_ref, emb_path_query does not exist')

        
    