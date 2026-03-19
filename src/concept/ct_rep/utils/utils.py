import anndata as ad
import scanpy as sc
import plotly.express as px
import umap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import issparse
from sklearn.decomposition import PCA
import re
from sklearn.model_selection import train_test_split

gene_name_to_id = pd.read_pickle("//p/scratch/cjinm16/dipippo1/scConcept/token_mapping/pc_gene_token_mapping.pkl")
gene_id_to_name = {v: k for k, v in gene_name_to_id.items()}

def is_raw_count(X):
    if isinstance(X, np.ndarray):
        return np.all(X >= 0) and (X == X.astype(int)).all()
    if issparse(X):
        return (X.data >= 0).all() and (X.data == X.data.astype(int)).all()
    
def normalize(X, method):
    EPSILON = 1e-3
    TOTAL_COUNT = 1e4
    
    assert isinstance(X, np.ndarray)
    if method == 'raw':
        X = X
    elif method == 'log1p':
        X = np.log1p(X)
    elif method == 'totalcount':
        X = X / (X.sum(axis=1)[:,None] + EPSILON) * TOTAL_COUNT
    elif method == 'totalcount/log1p' or method == 'totalcount_log1p':
        X = X / (X.sum(axis=1)[:,None] + EPSILON) * TOTAL_COUNT
        X = np.log1p(X)
    else:
        raise ValueError(f'Invalid normalization: {method}')
    return X


def extract_baseline_embs(adata_ref, adata_query, method, key = 'X_emb'):
    assert is_raw_count(adata_ref.X), 'adata_ref is not raw count'
    assert is_raw_count(adata_query.X), 'adata_query is not raw count'
    
    normalization = method.lower().split()[-1]
    print(f'Normalization: {normalization}')
    
    if 'inner' in method.lower():
        shared_genes = adata_ref.var.index.intersection(adata_query.var.index)
    elif 'outer' in method.lower():
        adata_query = adapt_vars(adata_query, adata_ref)
        shared_genes = adata_ref.var.index.intersection(adata_query.var.index)
    else:
        raise ValueError(f'Invalid method: {method}')

    X_ref = normalize(adata_ref[:, shared_genes].X.toarray(), normalization)
    X_query = normalize(adata_query[:, shared_genes].X.toarray(), normalization)

    if "raw" in method.lower():
        adata_ref.obsm[key] = X_ref
        adata_query.obsm[key] = X_query
        
    elif "pca" in method.lower():
        n_comp = int(re.search(r'\d+', method).group())
        if n_comp > X_ref.shape[1]:
            raise ValueError(f'n_comp ({n_comp}) > n_genes ({X_ref.shape[1]})')
        pca = PCA(n_components=n_comp)
        pca.fit(X_ref)
        adata_ref.obsm[key] = pca.transform(X_ref)
        adata_query.obsm[key] = pca.transform(X_query)
    return adata_ref, adata_query


def split_train_test(adata_train, adata_test, split_key_ref, split_key_query, split_value_ref, split_value_query):
    if split_value_ref=="%":
        assert (adata_train.obs.index == adata_test.obs.index).all(), 'index mismatch'
        train_split, test_split = train_test_split(list(adata_train.obs[split_key_ref].unique()), test_size=0.5, random_state=42)
        print(f'train split: {train_split}, test split: {test_split}')
        adata_train = adata_train[adata_train.obs[split_key_ref].isin(train_split)]
        adata_test = adata_test[adata_test.obs[split_key_query].isin(test_split)]
    elif split_value_ref=="*":
        adata_train = adata_train
        adata_test = adata_test
    else:
        adata_train = adata_train[adata_train.obs[split_key_ref].str.contains(split_value_ref)]
        adata_test = adata_test[adata_test.obs[split_key_query].str.contains(split_value_query)]
    return adata_train, adata_test


def get_discrete_palette(seed=42):
    # palette = sns.color_palette("colorblind", 20) + sns.color_palette("tab20", 20)
    palette = sns.color_palette("tab20", 20)
    palette = [palette[i] for i in range(0, 20, 2)] + [palette[i] for i in range(1, 20, 2)]
    palette = palette + sns.color_palette("tab20b", 20) # + sns.color_palette("tab20c", 20)
    rng = np.random.default_rng(seed=seed)
    rng.shuffle(palette)
    return palette

def reverse_normalization(X):
    
    assert isinstance(X, np.ndarray), 'X must be a numpy array'
    
    # Reverse the log1p normalization of the data
    X = np.expm1(X)
    
    # Replace zeros with np.nan to ignore them in the min calculation
    X_nonzero = np.where(X == 0, np.nan, X)

    # Find the minimum value for each row, ignoring np.nan
    min_nonzero_values = np.nanmin(X_nonzero, axis=1)

    # Assume that the minimum nonzero value for each cell is 1
    X = X / min_nonzero_values[:, None]

    # check that the result is close to an integer
    assert (X - np.round(X)).max() < 0.01, 'X is not close to an integer'

    return np.round(X).astype(int)


def get_hvgs(adata, n_top_genes=1000, flavor='seurat'):
    adata = adata.copy()
    if flavor not in ['seurat_v3', 'seurat_v3_paper']:
        print(f'using log data for flavor {flavor}')
        adata.layers['X_log1p'] = sc.pp.log1p(adata, copy=True).X
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor, layer='X_log1p')
    else:
        print(f'using raw data for flavor {flavor}')
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)
    return adata.var[adata.var.highly_variable].index.astype(str)

def adapt_vars(adata_target, adata_ref, return_mask=False):
    var_names_original = adata_target.var.index
    adata = ad.concat([adata_ref, adata_target], join='outer')
    adata_target = adata[len(adata_ref):]
    adata_target = adata_target[:, adata_ref.var.index]
    assert (adata_ref.var_names == adata_target.var_names).all()
    if not return_mask:
        return adata_target
    else:
        mask = adata_target.var.index.isin(var_names_original)
        return adata_target, mask
    

def sample_data(adata, key, n_obs=100):
    adatas = []
    values = adata.obs[key].unique()
    
    for value in values:
        adata_subsampled = adata[adata.obs[key] == value]
        if len(adata_subsampled) > n_obs:
            adata_subsampled = sc.pp.subsample(adata_subsampled, n_obs=n_obs, copy=True)
        adatas.append(adata_subsampled)
    
    adata = ad.concat(adatas, join='outer', merge='same')
    
    return adata



def get_panel(panel_name):
    panel_dir = Path('/lustre/groups/ml01/projects/contrastive_transformer/panels/human')
    filename = panel_name + '.csv'
    panel = pd.read_csv(panel_dir / filename)
    panel = panel['Ensembl_ID'].rename('feature_id', inplace=True).values
    return panel

def plot_umap(adata, color_keys, use_rep='X_emb', max_colors: int | None = 10, metric='euclidean', subsample=5000, n_neighbors=15, train_obs_key=None, train_value=None, color_discrete_map=None, title=None, fig_save_path=None, fig_save_name=None):
    
    if 'X_umap' in adata.obsm:
        del adata.obsm['X_umap']
    if 'umap' in adata.obsm:
        del adata.obsm['umap']
    
    for color_key in color_keys:
        if max_colors is not None and adata.obs[color_key].dtype == 'category':
            print(f'taking top {max_colors} {color_key}')
            adata = adata[adata.obs[color_key].isin(adata.obs[color_key].value_counts().index[:max_colors])]
    
    if subsample is not None:
        adata = sc.pp.subsample(adata, n_obs=min(subsample, len(adata)), copy=True)
    
    if train_obs_key is not None and train_value is not None:
        adata_train = adata[adata.obs[train_obs_key] == train_value]
        umap_model = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.5, spread=1.0, metric=metric, random_state=0).fit(adata_train.obsm[use_rep])
        adata.obsm['X_umap'] = umap_model.transform(adata.obsm[use_rep])
    else:
        sc.pp.neighbors(adata, use_rep=use_rep, metric=metric, n_neighbors=n_neighbors)
        sc.tl.umap(adata)
        
    for color_key in color_keys:
        fig = get_umap(adata, color_key, title, color_discrete_map)
        if fig_save_path is not None and fig_save_name is not None:
            fig.write_image(Path(fig_save_path) / f'{fig_save_name}_{color_key}.pdf')
        fig.show()


def get_umap(adata, color_key, title, color_discrete_map, marker_size=3, marker_opacity=1.0, only_plot=False):
    fig = px.scatter(x=adata.obsm['X_umap'][:,0], y=adata.obsm['X_umap'][:,1], color=adata.obs[color_key], title=title,
                        color_discrete_map=color_discrete_map, size_max=1,
                        color_continuous_scale=px.colors.sequential.Turbo)
    fig.update_layout(
        width=700 , #+ (not only_plot) * 8.4 * adata.obs[color_key].astype(str).str.len().max(), 
        height=700,
        legend=dict(itemsizing='constant', font=dict(size=14)),  # Increase legend marker size
        showlegend=not only_plot,  # Hide legend
        xaxis=dict(showticklabels=False, title=''),  # Remove x-axis ticks and title
        yaxis=dict(showticklabels=False, title=''),  # Remove y-axis ticks and title
    )
    fig.update_traces(marker_size=marker_size, marker_opacity=marker_opacity)
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(
        mirror=True,
        ticks=None,
        showline=True,
        linecolor='black',
        gridcolor='white'
    )
    fig.update_yaxes(
        mirror=True,
        ticks=None,
        showline=True,
        linecolor='black',
        gridcolor='white',
        scaleanchor = "x",
        scaleratio = 1
    )
    return fig

def get_umap_scanpy(adata, color_key, title=None, color_discrete_map='viridis', marker_size=3, marker_opacity=1.0, only_plot=False):
    """
    Plot UMAP using Scanpy's sc.pl.umap.

    Parameters:
        adata: AnnData object
        color_key: str - column in adata.obs to color by
        title: str - plot title
        color_map: str or matplotlib colormap - colormap for continuous variables
        marker_size: float - size of points
        marker_opacity: float - alpha transparency of points (not directly supported by sc.pl.umap)
        only_plot: bool - if True, hides legend and title
    """
    fig = sc.pl.umap(
        adata,
        color=color_key,
        title='' if only_plot else title,
        size=marker_size,
        show=False,
        palette=color_discrete_map,
        legend_loc=None if only_plot else 'right margin',
        frameon=True,
        return_fig=True,
    )
    return fig

def compare_feature_distribution(adata_1, adata_2, shared_features=False, min_count=1):
    from ct_rep.utils.utils import normalize
    if shared_features:
        shared_features = adata_1.var.index.intersection(adata_2.var.index)
        adata_1 = adata_1[:, shared_features]
        adata_2 = adata_2[:, shared_features]
    
    n_obs = min(len(adata_1), len(adata_2), 5000)
    adata_1_plot = sc.pp.subsample(adata_1, n_obs=n_obs, copy=True)
    adata_2_plot = sc.pp.subsample(adata_2, n_obs=n_obs, copy=True)

    adata_1_plot.X = normalize(adata_1_plot.X.toarray(), 'raw')
    adata_2_plot.X = normalize(adata_2_plot.X.toarray(), 'raw')

    # adata_1_plot.X = np.apply_along_axis(normalize, 1, adata_1_plot.X, 'binning/20')
    # adata_2_plot.X = np.apply_along_axis(normalize, 1, adata_2_plot.X, 'binning/20')

    # min_count = 1
    # bins = np.linspace(min_count, 50, 50)
    # pd.Series(adata_2_plot.X.flatten()).hist(bins=bins, color='orange', alpha=0.5)
    # plt.ylabel('Frequency')
    # plt.title(f'counts {dataset_test}')

    plt.figure()
    data_1 = pd.Series((adata_1_plot.X>0).sum(axis=1).squeeze())
    data_2 = pd.Series((adata_2_plot.X>0).sum(axis=1).squeeze())
    bins = np.linspace(0, max(data_1.max(), data_2.max()) , 50)
    data_1.hist(bins=bins, color='blue', alpha=0.5, label='adata_1')
    data_2.hist(bins=bins, color='red', alpha=0.5, label='adata_2')
    plt.legend()
    plt.title(f'adata_1 max value:{data_1.max()}, adata_2 max value:{data_2.max()}')
    plt.ylabel('Frequency')
    plt.xlabel(f'Number of nonzero features')
    
    plt.figure()
    data_1 = pd.Series((adata_1_plot.X>0).mean(axis=1).squeeze())
    data_2 = pd.Series((adata_2_plot.X>0).mean(axis=1).squeeze())
    bins = np.linspace(0, 1, 20)
    data_1.hist(bins=bins, color='blue', alpha=0.5, label='adata_1')
    data_2.hist(bins=bins, color='red', alpha=0.5, label='adata_2')
    plt.legend()
    plt.title(f'adata_1 max value:{data_1.max()}, adata_2 max value:{data_2.max()}')
    plt.ylabel('Frequency')
    plt.xlabel(f'Ratio of nonzero features')

    plt.figure()
    data_1 = pd.Series(adata_1_plot.X.sum(axis=1).squeeze())
    data_2 = pd.Series(adata_2_plot.X.sum(axis=1).squeeze())
    bins = np.linspace(0, max(data_1.max(), data_2.max()) , 50)
    data_1.hist(bins=bins, color='blue', alpha=0.5, label='adata_1')
    data_2.hist(bins=bins, color='red', alpha=0.5, label='adata_2')
    plt.legend()
    plt.title(f'adata_1 max value:{data_1.max()}, adata_2 max value:{data_2.max()}')
    plt.ylabel('Frequency')
    plt.xlabel(f'Sum(counts)')

    plt.figure()
    data_1 = adata_1_plot.X.flatten()
    data_2 = adata_2_plot.X.flatten()
    data_1 = data_1[data_1>=min_count]
    data_2 = data_2[data_2>=min_count]
    n_nz = min(len(data_1), len(data_2), 5000)
    data_1 = pd.Series(np.random.choice(data_1, n_nz))
    data_2 = pd.Series(np.random.choice(data_2, n_nz))
    bins = np.linspace(min_count, max(data_1.max(), data_2.max()) , 50)
    data_1.hist(bins=bins, color='blue', alpha=0.5, label='adata_1')
    data_2.hist(bins=bins, color='red', alpha=0.5, label='adata_2')
    plt.legend()
    plt.ylabel('Frequency')
    plt.title(f'adata_1 max value:{data_1.max()}, adata_2 max value:{data_2.max()}')
    plt.xlabel(f'Nonezero values')


def matching_acc(adata_1, adata_2, batch_size=256, union=True):
    batch_sizes = []
    acc_conts = []
    top5_acc_conts = []
    for i in range(300):
        random_indices = np.random.choice(len(adata_1), size=batch_size, replace=False)
        
        if union:
            cell_embs_11 = torch.tensor(adata_1[random_indices[:batch_size//2]].obsm['X_emb'])
            cell_embs_12 = torch.tensor(adata_2[random_indices[batch_size//2:]].obsm['X_emb'])
            cell_embs_1 = torch.cat([cell_embs_11, cell_embs_12])
            cell_embs_21 = torch.tensor(adata_2[random_indices[:batch_size//2]].obsm['X_emb'])
            cell_embs_22 = torch.tensor(adata_1[random_indices[batch_size//2:]].obsm['X_emb'])
            cell_embs_2 = torch.cat([cell_embs_21, cell_embs_22])
        else:
            cell_embs_1 = torch.tensor(adata_1[random_indices[:batch_size]].obsm['X_emb'])
            cell_embs_2 = torch.tensor(adata_2[random_indices[:batch_size]].obsm['X_emb'])
        
        logits = torch.mm(cell_embs_1, cell_embs_2.t())
        cont_target = torch.arange(len(logits))
        acc_cont = (logits.argmax(dim=1) == cont_target).float().mean()
        top5_acc_cont = (logits.topk(5, dim=1)[1] == cont_target.unsqueeze(1)).any(dim=1).float().mean()
        batch_sizes.append(len(logits))
        acc_conts.append(acc_cont)
        top5_acc_conts.append(top5_acc_cont)
    print(f'union:{union}, batch_size:{np.mean(batch_sizes)}, acc_cont: {np.mean(acc_conts)*100:.3f}, top5_acc_cont: {np.mean(top5_acc_conts)*100:.3f}')
    return np.mean(acc_conts)


def same_batch_ratio(adata, batch, emb_key='X_emb', n_neighbors=100, metric='cosine', query_key=None, query_value=None):
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric=metric).fit(adata.obsm[emb_key])

    # Select a random sample from adata
    if query_key is not None:
        indices = np.where(adata.obs[query_key] == query_value)[0]
        print('indices: ', len(indices))
    else:
        indices = np.arange(len(adata))
        
    size = min(len(indices), 10000)
    random_index = np.random.choice(indices, size=size, replace=False)
    random_sample = adata.obsm[emb_key][random_index]
    random_batch = adata.obs.iloc[random_index][batch].values.reshape(-1, 1)

    distances, indices = nbrs.kneighbors(random_sample)

    # Count how many of the closest neighbors are from spatial
    neighbor_batches = adata.obs.iloc[indices.flatten()][batch].values.reshape(indices.shape)
    same_batch_ratio = (neighbor_batches == random_batch).mean(axis=1).mean()

    print(f'Same {batch} Ratio: {same_batch_ratio*100:.3f}')
    return same_batch_ratio

def shared_feature_stats(adata, batch_size=256):
    random_indices = np.random.choice(len(adata), size=batch_size, replace=False)
    num_shared_featrues = []
    for i in random_indices:
        for j in random_indices:
            if i < j :
                num_shared_featrues.append(np.intersect1d(adata.X[i].nonzero()[1], adata.X[j].nonzero()[1]).size / np.union1d(adata.X[i].nonzero()[1], adata.X[j].nonzero()[1]).size)
    
    print(f'Average % of shared features: %{np.median(num_shared_featrues)*100:.3f}')
    return np.array(num_shared_featrues)



def get_unbalanced_cm(true_labels, y_pred, classes_order=None):
    # Compute confusion matrix manually, allowing for different sets of true and predicted classes
    true_classes_ = np.unique(true_labels)
    pred_classes_ = np.unique(y_pred)
    if classes_order is not None:
        true_classes_ = [c for c in classes_order if c in true_classes_]
        pred_classes_ = [c for c in classes_order if c in pred_classes_]
        
    cm = np.zeros((len(true_classes_), len(pred_classes_)), dtype=int)
    true_label_to_index = {label: idx for idx, label in enumerate(true_classes_)}
    pred_label_to_index = {label: idx for idx, label in enumerate(pred_classes_)}
    for true_label, pred_label in zip(true_labels, y_pred):
        if true_label in true_label_to_index and pred_label in pred_label_to_index:
            i = true_label_to_index[true_label]
            j = pred_label_to_index[pred_label]
            cm[i, j] += 1
    return cm, true_classes_, pred_classes_
    
