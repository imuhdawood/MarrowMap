import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import seaborn as sns
import os

import sys
sys.path.append('.')

from utils.utils_io import mkdir
from configuration import STConfig

# Define constants and paths
K_RANGE = list(range(2, 6))
K = 3
FIG_SIZE = (8, 10)

cfg = STConfig()
PATH = f'{cfg.pth_downstream_out}/obj_w_Cell_annotations_RC_22122024_final_with_osteo/BandFeatures'
OUTPATH = f'{cfg.pth_downstream_out}/obj_w_Cell_annotations_RC_22122024_final_with_osteo/BandFeatureCluster'

# Getting optimal number of clusters using the Silhouette method
def determine_optimal_clusters(data, features_col):
    silhouette_scores = []
    for k in K_RANGE:
        kmeans = KMeans(n_clusters=k, random_state=42)
        labels = kmeans.fit_predict(data[features_col])
        silhouette_avg = silhouette_score(data[features_col], labels)
        silhouette_scores.append(silhouette_avg)
    
    # Plot the Silhouette method graph
    plt.figure(figsize=FIG_SIZE)
    plt.plot(K_RANGE, silhouette_scores, 'bo-')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Method For Optimal K')
    plt.savefig(f'{PLOTS_DIR}/silhouette_method.png')

    optimal_k = K_RANGE[np.argmax(silhouette_scores)]
    return optimal_k

# Running KMeans clustering to identify spatiotypes of each cell type
def perform_clustering_and_generate_plots(data, features_col, cell_type,celltypes_selected=None):
    optimal_k = K#determine_optimal_clusters(data, features_col)
    if cell_type in celltypes_selected:
        optimal_k = celltypes_selected[cell_type]
    kmeans = KMeans(n_clusters=optimal_k, random_state=42)
    data['cluster'] = kmeans.fit_predict(data[features_col])
    
    cluster_means = data[features_col + ['cluster']].groupby('cluster').mean() + 1e-6
    reorder = sorted(cluster_means.columns)
    cluster_means = cluster_means.loc[:,reorder]

    # Spatial cellular abundnace across spatiotypes
    plt.figure(figsize=FIG_SIZE)
    sns.heatmap(cluster_means.T, cmap="YlGnBu", annot=True, vmin=0, vmax=1)
    plt.title('Mean Microenvironmental Features per Spatiotype')
    plt.xlabel('Cluster')
    plt.ylabel('Feature')
    plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/{cell_type}_abundance.png')

    # Synthetic reference cluster and log2 fold change calculation
    synthetic_reference = pd.Series(1 / len(cluster_means.columns), index=cluster_means.columns)
    cluster_means.loc['synthetic_reference'] = synthetic_reference
    log2_fold_change = np.log2(cluster_means / cluster_means.loc['synthetic_reference'])
    log2_fold_change.replace([np.inf, -np.inf], np.nan, inplace=True)
    log2_fold_change = log2_fold_change.drop('synthetic_reference')

    # Visualize the log2 fold change using a heatmap
    plt.figure(figsize=FIG_SIZE)
    sns.heatmap(log2_fold_change.T, cmap="RdBu_r", annot=True, center=0)
    plt.title('Log2 Fold Change of Microenvironmental Features per Cluster (Relative to Synthetic Reference)')
    plt.xlabel('Cluster')
    plt.ylabel('Feature')
    plt.tight_layout()
    plt.savefig(f'{PLOTS_DIR}/{cell_type}_abundance_log_fold.png')

    # Save clustered data to CSV
    data.to_csv(f'{CSV_DIR}/{cell_type}.csv')

# Main script execution
if __name__ == "__main__":

    celltypes_selected = {
        'Stromal': 3, 'Lymphocyte': 3, 'HSPC': 3,
        'Endothelial': 3, 'Immature_myeloid': 3, 'Myeloid': 3,
        'Erythroid': 3, 'MNP': 3,
        'Megakaryocyte': 3, 'Granulocyte_mast': 3, 'Osteo-MSC':3
    }

    for cell_type in celltypes_selected:
        if cell_type in ['Granulocyte/mast']:
            cell_type = 'Granulocyte_mast'
        INPUT_FILE = f'{PATH}/band_features_{cell_type}_cell_id_br2.csv'
        OUT_PATH = f'{OUTPATH}/br2_k{K}/{cell_type}'
        #OUT_PATH = f'{OUTPATH}/br2_custom/{cell_type}'
        mkdir(OUT_PATH)
        PLOTS_DIR = f'{OUT_PATH}/plots'
        CSV_DIR = f'{OUT_PATH}/csv'
        mkdir(PLOTS_DIR)
        mkdir(CSV_DIR)
        print(f'Processing cell type {cell_type}')
        data = pd.read_csv(INPUT_FILE)
        # idx = np.random.randint(0,len(data),500)
        # data = data.iloc[idx,:]
        features_col = list(set(data.columns) - {'sample_id', 'label','cell_id'})
        perform_clustering_and_generate_plots(data, features_col, cell_type,celltypes_selected = celltypes_selected)
