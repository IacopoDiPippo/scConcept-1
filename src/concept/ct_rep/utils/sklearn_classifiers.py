import argparse
import numpy as np
import scanpy as sc
from sklearn.preprocessing import LabelEncoder
from sklearn.linear_model import SGDClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, f1_score, mean_squared_error, root_mean_squared_error
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from scipy import stats

def pearson_corr_score(y_true, y_pred, normalize=False):
    if normalize:
        y_true = y_true / (y_true.sum(axis=1)[:, None] + 0.001)
        y_pred = y_pred / (y_pred.sum(axis=1)[:, None] + 0.001)
    if y_true.shape != y_pred.shape: raise ValueError("Shapes are different.")
    corrsum = 0
    corrs = []
    for i in range(y_true.shape[1]):
        corr = stats.pearsonr(y_true[:, i], y_pred[:, i])[0]
        corrsum += corr if not np.isnan(corr) else 0
        corrs.append(corr if not np.isnan(corr) else 0)
    return corrsum / y_true.shape[1], np.array(corrs)


def spearman_corr_score(y_true, y_pred):    
    if y_true.shape != y_pred.shape: raise ValueError("Shapes are different.")
    corrsum = 0
    for i in range(y_true.shape[1]):
        corr = stats.spearmanr(y_true[:, i], y_pred[:, i])[0]
        corrsum += corr if not np.isnan(corr) else 0
    return corrsum / y_true.shape[1]


def knn_on_embeddings(train_embs, test_embs, train_labels, test_labels, n_neighbors=10, metric='minkowski', pred_labels=None):
    print(f'Fitting KNN ({metric}) with n_neighbors: ', n_neighbors)
    knn = KNeighborsClassifier(n_neighbors=n_neighbors, metric=metric, algorithm='brute')
    if pred_labels is not None:
        mask = np.isin(train_labels, pred_labels)
        train_labels = train_labels[mask]
        train_embs = train_embs[mask]
    knn.fit(train_embs, train_labels)
    y_pred = knn.predict(test_embs)
    return y_pred

def logistic_regression_on_embeddings(train_embs, test_embs, train_labels, test_labels, max_iter=1000):
    logistic_regression = SGDClassifier(loss='log_loss', max_iter=max_iter, n_jobs=5, verbose=True)
    logistic_regression.fit(train_embs, train_labels)
    y_pred = logistic_regression.predict(test_embs)
    return y_pred


def mlp_on_embeddings(train_embs, test_embs, train_labels, test_labels, max_iter=200, hidden_layer_sizes=[128, 128], early_stopping=False):    
    clf = MLPClassifier(max_iter=max_iter, hidden_layer_sizes=hidden_layer_sizes, verbose=True, early_stopping=early_stopping, tol=0.001, n_iter_no_change=5, learning_rate_init=0.0001)
    clf.fit(train_embs, train_labels)
    y_pred = clf.predict(test_embs)
    return y_pred

def knn_reg_on_embeddings(train_embs, test_embs, train_targets, n_neighbors=10, metric='minkowski'):
    print(f'Fitting KNN ({metric}) with n_neighbors: ', n_neighbors)
    knn = KNeighborsRegressor(n_neighbors=n_neighbors, metric=metric, algorithm='brute')
    knn.fit(train_embs, train_targets)
    y_pred = knn.predict(test_embs)
    return y_pred
