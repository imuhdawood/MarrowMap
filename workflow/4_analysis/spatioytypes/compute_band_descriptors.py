import json
import os
import multiprocessing as mp
from tqdm import tqdm
import pandas as pd
import numpy as np
import anndata as ad

from sklearn.neighbors import KDTree

import sys
sys.path.append('.')

from utils.utils_io import mkdir
from configuration import STConfig


def compute_microenvironment(cell_idx, radii, data, tree, cell_bands_template, mapping):
    center = data[['x', 'y']].loc[cell_idx].values.reshape(1, -1)
    
    # Initialize a dictionary to hold counts
    cell_bands = cell_bands_template.copy()

    previous_radius = 0
    for r in radii:
        # Query cells within the current radius
        idx_within_radius = tree.query_radius(center, r=r)[0]

        # Get only the cells in the current band, not overlapping with previous radius
        if previous_radius > 0:
            idx_within_previous = tree.query_radius(center, r=previous_radius)[0]
            idx_in_band = np.setdiff1d(idx_within_radius, idx_within_previous)
        else:
            idx_in_band = idx_within_radius

        if len(idx_in_band) > 0:
            # Count the occurrences of each cell type within the band
            idx_in_band = [mapping[id] for id in idx_in_band]
            type_counts = data.loc[idx_in_band]['cell_type'].value_counts(normalize=True)
            for cell_type, count in type_counts.items():
                cell_bands[f'{cell_type}_{r}'] = count
        previous_radius = r
    cell_bands['cell_id'] = data.loc[cell_idx].cell_id    
    return cell_bands

def process_sample(sample_id, consolDf, cell_of_interest, celltypes_selected, radii):

    fDf = consolDf[consolDf.sample_key.isin([sample_id])]
    data = fDf[['x', 'y','cell_type','cell_id']]
    data = data[data.cell_type.isin(celltypes_selected)]
    cell_bands_template = {f'{c}_{r}': 0 for c in celltypes_selected for r in radii}
    # Build KD-Tree for spatial queries
    tree = KDTree(data[['x', 'y']].values)
    # Define the safe boundary based on the maximum radius
    max_radius = max(radii)
    x_min_safe, y_min_safe = max_radius, max_radius
    x_max_safe, y_max_safe = data['x'].max() - max_radius, data['y'].max() - max_radius
    mapping = dict(
        zip(
            list(range(data.shape[0])),
            data.index.tolist()
    ))
    # Filter only the cells of interest within the safe boundary
    no_filtering = data['cell_type'] == cell_of_interest
    print('Without boundary check',sum(no_filtering))
    subset_idx = data[
        (data['cell_type'] == cell_of_interest) &
        (data['x'] > x_min_safe) & (data['x'] < x_max_safe) &
        (data['y'] > y_min_safe) & (data['y'] < y_max_safe)
    ].index.to_numpy()
    print('After filtering', len(subset_idx))
    
    # Parallelize the microenvironment computation
    with mp.Pool(mp.cpu_count()) as pool:#
        features = pool.starmap(
            compute_microenvironment, 
            [(idx, radii, data, tree, cell_bands_template, mapping) for idx in subset_idx]
        )
    
    features_df = pd.DataFrame(features)
    features_df['sample_id'] = sample_id

    return features_df

def main():

    skey = 'sample_key'
    cell_type_col = 'obj.anno_5_w_megs_w_stromal'

    cfg = STConfig()
    meta_df = pd.read_csv(cfg.pth_meta_csv)
    samples = meta_df[skey].unique().tolist()
    print(f'Number of samples {len(samples)}')

    # reading consolidated adata object
    adata = ad.read_h5ad(cfg.pth_consol_adata)
    adata = adata[adata.obs["cell_status"] != "adipocyte"]

    # reading cell annotation 
    annot_df = pd.read_csv(cfg.pth_cell_annotations_final)
    annot_df[skey] = annot_df['cell_id'].str[-8:]
    annot_df = annot_df.loc[:, [skey,'cell_id' ,cell_type_col]]
    annot_df.index = annot_df['cell_id'].tolist()
    annot_df.drop(columns=['cell_id'], inplace=True)

    adata.obs = adata.obs.merge(annot_df, left_index = True, right_index = True)
    adata.obs[cell_type_col] = adata.obs[cell_type_col].replace({'Adipo-MSC':'Stromal', 'SMC':'Stromal'})
    mask = (adata.obs['annotation'] == 'keep') | (adata.obs[cell_type_col].notna())
    adata = adata[mask]
    consolDf = adata.obs.copy()
    consolDf[['x', 'y']] = np.array(adata.obsm['spatial'])
    consolDf.rename(columns = {cell_type_col:'cell_type'}, inplace=True)
    consolDf = consolDf[~consolDf.cell_type.isin(['remove','CD69'])]

    selected_columns = ['cell_type','x','y',skey,'cell_id']

    OUTPUT_DIR = f'{cfg.pth_downstream_out}/BandFeatures'
    mkdir(OUTPUT_DIR)

    celltypes_selected = consolDf['cell_type'].unique()
    radii = [100, 200, 300, 400, 500]# These are band radii 1 clusters
    for cell_of_interest in celltypes_selected:
        print('Processing *******************************************')
        print('********************************************************')
        print(cell_of_interest)
        allSampleDf = pd.DataFrame()

        for sample_id in tqdm(samples):
            mDf = consolDf[consolDf[skey].isin([sample_id])].loc[:,selected_columns]
            print(mDf.dropna().shape, sample_id)

            features_df = process_sample(sample_id,consolDf, cell_of_interest, celltypes_selected, radii)
            allSampleDf = pd.concat([allSampleDf, features_df], ignore_index=True)
        if cell_of_interest in ['Granulocyte/mast']:
            cell_of_interest = 'Granulocyte_mast'
        output_file = f'{OUTPUT_DIR}/band_features_{cell_of_interest}_cell_id_br2.csv'
        allSampleDf.to_csv(output_file, index=False)
        print(f'Saved results to {output_file}')

if __name__ == "__main__":
    main()
