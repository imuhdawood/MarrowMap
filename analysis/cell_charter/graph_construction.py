import os
from tqdm import tqdm
import spatialdata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import anndata as ad
import squidpy as sq
import networkx as nx
from geopandas import GeoDataFrame
from scipy.sparse import csr_matrix

import sys
sys.path.append('..')
from configuration import STConfig
from wrappers.sdata_initalizer import SDataInitalizer
from wrappers.sdata_customizer import SDataCustomizer

from utils.utils_colors import cn_color_palette
from utils.utils_geom import align_sample_polygons

cfg = STConfig()
adata = ad.read_h5ad(cfg.pth_consol_adata)

cfg = STConfig()
meta_df = pd.read_csv(cfg.pth_meta_csv)
samples = meta_df['sample_key'].tolist()

for sample_id in samples:
    customizer = SDataCustomizer(cfg, sample_key=sample_id)
    sdata_file_path = f'{cfg.pth_sdata}/{sample_id}_no_he.zarr'
    sdata_obj = spatialdata.read_zarr(sdata_file_path)
    # selecting cells falling in IT regions of current sample
    adata_sub = adata[adata.obs['cell_id'].str.contains(sample_id)].copy()
    adata_sub_itr = adata_sub[adata_sub.obs['it_regions'] != "non_intertrabecular"]
    mask = ~np.isnan(adata_sub_itr.obsm["spatial"]).any(axis=1)
    adata_sub_itr = adata_sub_itr[mask].copy()
    cell_status_df = adata_sub_itr.obs[['extra_new_cell_status', 'meg_phenotype']]

    # Selecting common cells in sdata_object and adata
    gdf_cells = sdata_obj["new_cell_boundaries"].copy()
    merged_gdf = gdf_cells.merge(cell_status_df, left_index=True, right_index=True)

    gdf_art = merged_gdf.loc[merged_gdf['extra_new_cell_status'] == 'artefact', ['meg_phenotype', 'geometry']].copy()
    gdf_meg = merged_gdf.loc[merged_gdf['meg_phenotype'] != "non_meg", ['meg_phenotype', 'geometry']].copy()

    # getting touching polygon artefact and non-megs to merge it for better graph
    gdf_art["geometry"] = gdf_art.buffer(10)
    gdf_meg["geometry"] = gdf_meg.buffer(10)
    sjoined = gdf_meg.sjoin(gdf_art, how="inner", predicate="intersects")

    meg_idices = sjoined.index.tolist()
    arf_idices = sjoined.index_right.tolist()

    # merge the touching artefact boundary to each meg
    for idx, meg_idx in enumerate(meg_idices):
        gdf_meg.loc[meg_idx, "geometry"] = gdf_meg.loc[meg_idx, "geometry"].union(gdf_art.loc[arf_idices[idx], "geometry"])

    # shrinking back the geometry to original bounds
    gdf_meg["geometry"] = gdf_meg.buffer(-10)

    # Getting a copy of cell geometry from sdata object 
    gdf_cells_new = sdata_obj["new_cell_boundaries"].copy()
    new_gdf = gdf_cells_new.drop(gdf_art.index)

    new_gdf.loc[gdf_meg.index, "geometry"] = gdf_meg["geometry"]

    # Selecting common index we are using for graph construction
    final_indexes = new_gdf.index.intersection(adata_sub_itr.obs.index.tolist())
    print(f"Total Cell Selected {len(final_indexes)}")
    new_gdf = new_gdf.loc[final_indexes]

    # Defining graph connectivity
    new_gdf = new_gdf.buffer(2)
    gdf = GeoDataFrame(new_gdf, columns=["geometry"], index=new_gdf.index)
    # self joint to identify touching polygon
    sjoined = gdf.sjoin(gdf, how="left", predicate="intersects")
    print(f"# of cell id {len(final_indexes)}, merged df shape {sjoined.shape}")

    # getting cell centroid as node coordiantes and connectivity from touching cells
    index_map = {key: idx for idx, key in enumerate(final_indexes)}
    edges = {
        (ix, iy) if ix < iy else (iy, ix)
        for x, y in zip(sjoined.index, sjoined.index_right)
        if x != y
        for ix, iy in [(index_map.get(x, -1), index_map.get(y, -1))]
    }
    centroids = gdf.geometry.centroid
    coords = np.column_stack((centroids.x.values, centroids.y.values))

    graph_adata = adata_sub_itr[final_indexes]

    graph_adata.uns["spatial_neighbors"] = {}
    graph_adata.uns["spatial_neighbors"]["params"] = {}

    # the steps below could be done even without nx.Graph but keeping it for debugging and easy visualization
    G = nx.Graph()
    pos_dict = {i: tuple(coords[i]) for i in range(coords.shape[0])}
    G.add_edges_from(edges)
    G.add_nodes_from(pos_dict)

    # In our case both wighted and unweighted adjacency matrix is same just for the sake of completion
    adj_matrix_unweighted = nx.to_scipy_sparse_array(G, weight=None)
    adj_matrix_unweighted = csr_matrix(adj_matrix_unweighted)

    adj_matrix_weighted = nx.to_scipy_sparse_array(G, weight='weight')
    adj_matrix_weighted = csr_matrix(adj_matrix_weighted)
    #add them back into the adata object

    adata.obsp['spatial_connectivities'] = adj_matrix_unweighted
    adata.obsp['spatial_distances'] = adj_matrix_weighted

    #saving adata with graph
    adata.write(os.path.join(cfg.pth_graph_adata, sample_id + "_sp.h5ad"))
