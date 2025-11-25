from sklearn.preprocessing import MinMaxScaler
from shapely import wkt
import os
from tiatoolbox.annotation.storage import SQLiteStore, Annotation
import numpy as np
import pandas as pd
from tqdm import tqdm

from tiatoolbox.wsicore.wsireader import OpenSlideWSIReader
# from tiatoolbox.tools.patchextraction import PatchExtractor

from typing import List
from shapely.geometry import Polygon
from sklearn.metrics import roc_auc_score


environ_feats = [
    'Endothelial-0', 'Endothelial-1', 'Endothelial-2',
    'Erythroid-0', 'Erythroid-1', 'Erythroid-2',
    'Granulocyte/mast-0', 'Granulocyte/mast-1', 'Granulocyte/mast-2',
    'HSPC-0', 'HSPC-1', 'HSPC-2',
    'Immature_myeloid-0', 'Immature_myeloid-1', 'Immature_myeloid-2',
    'Lymphocyte-0', 'Lymphocyte-1', 'Lymphocyte-2',
    'MNP-0', 'MNP-1', 'MNP-2',
    'Megakaryocyte-0', 'Megakaryocyte-1', 'Megakaryocyte-2',
    'Myeloid-0', 'Myeloid-1', 'Myeloid-2',
    'Stromal-0', 'Stromal-1', 'Stromal-2',
    'Osteo-MSC-0', 'Osteo-MSC-1', 'Osteo-MSC-2'
    ]

# Selected set of features based on statistical significance
environ_feats_selected = [
     'Erythroid-1', 'Erythroid-2',                                                           
    'Immature_myeloid-0', 'Immature_myeloid-2',
     'MNP-2',
    'Megakaryocyte-0', 'Megakaryocyte-1',
    'Myeloid-0',
    'Stromal-1',
    'Osteo-MSC-2'
    ]

# celltypes_feats = [
#     'Immature_myeloid', 'Lymphocyte', 'Stromal',
#     'MNP', 'Endothelial', 'Erythroid',
#     'Megakaryocyte', 'HSPC', 'Granulocyte/mast'
#                     ]


celltypes_feats = [
    'Stromal', 'Lymphocyte', 'HSPC',
    'Endothelial', 'Immature_myeloid', 'Myeloid',
    'Erythroid', 'MNP', 'Osteo-MSC',
    'Megakaryocyte','Granulocyte/mast'
]

    # Annot 3
celltypes_selected_anno3 = [
    'DC', 'Erythroid', 'Myeloid',
    'Monocyte', 'B_cell', 'T_cell',
    'Megakaryocyte', 'Endothelial', 'Osteo-MSC',
    'HSPC', 'Macrophage', 'Plasma_cell',
    'GMP', 'SMC', 'Stromal',
    'Adipo-MSC', 'Granulocyte/mast'
]

# Annot 1
celltypes_selected_anno1 = [
    'DC', 'Erythroid', 'Myeloid', 'Non-classical_monocyte', 'B_cell',
    'CD8_T_cell', 'Megakaryocyte', 'Endothelial', 'IL7R_expressing',
    'Osteo-MSC', 'HSPC', 'Macrophage',
    'Pro-inflammatory_classical_monocyte', 'Plasma_cell', 'GMP',
    'Cycling_immature_erythroid', 'Cytotoxic_T_cell/NK', 'SMC',
    'Stromal', 'Cycling_myeloid', 'Adipo-MSC', 'Classical_monocyte',
    'Granulocyte/mast', 'CD69_expressing', 'T_cell',
#    'Granylocyte/mast'
    ]

def min_max_scaler(vector):
    min_val = vector.min()
    max_val = vector.max()

    # Perform min-max scaling to range [0, 1]
    scaled_vector = (vector - min_val) / (max_val - min_val)
    return scaled_vector

def poly_to_area(poly_str):
    polygon = wkt.loads(poly_str)
    area = polygon.area
    return area

def mkdirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def write_annotations_its_level(
        itr_df,
        feats_cols, 
        OUT_ANNOT_DIR,
        scaling=True, 
        specific_sample_ids=None
    ):

    os.makedirs(OUT_ANNOT_DIR, exist_ok=True)

    if scaling and 'scores' in feats_cols:
        itr_df[['scores']] = MinMaxScaler().fit_transform(itr_df[['scores']])
        itr_df[itr_df[feats_cols]>1] = 1
    
    # Create unique set of sample IDs from the index
    samples_set = list(set(itr_df.index.str.split('_').str[:2].str.join('_')))
    
    for sample_id in samples_set:
        subset = itr_df[itr_df.index.str.contains(sample_id)]
        annotations = []
        
        # Optional specific sample ID handling for Error Handling
        if specific_sample_ids and sample_id in specific_sample_ids:
            import pdb ;pdb.set_trace()
            print(f"Processing specific sample ID: {sample_id}")
        
        SQ = SQLiteStore()
        print('Writing annotation file')
        
        # Generate annotations for each row in the subset
        for itdx in subset.index:
            props = subset.loc[itdx, feats_cols].astype('float').to_dict()
            geometry = subset.loc[itdx, 'geometry']
            polygon = wkt.loads(geometry)
            
            annotations.append(
                Annotation(
                    geometry=polygon,
                    properties=props
                )
            )
        
        # Save annotations to database
        SQ.append_many(annotations)
        out_db_file = f'{OUT_ANNOT_DIR}/{sample_id}.db'
        SQ.create_index("area", '"area"')
        SQ.dump(out_db_file)


def extract_patches_from_wsi(
    cases, 
    WSI_DIR,
    resolution = 0.2125, 
    patch_shape=(1024, 1024), 
    stride_shape=(512, 512)
    ):
    """
    Extracts patches from whole-slide images (WSIs) and returns patch coordinates in a DataFrame.
    
    Args:
        cases (list): List of sample IDs to process.
        WSI_DIR (str): Directory where WSI files are stored.
        PatchExtractor (class): PatchExtractor class with a get_coordinates method.
        OpenSlideWSIReader (class): OpenSlideWSIReader class to read WSIs.
        resolution (float): Image resolution, default is 0.2125.
        patch_shape (tuple): Shape of each patch, default is (1024, 1024).
        stride_shape (tuple): Stride shape for patch extraction, default is (512, 512).
    
    Returns:
        pd.DataFrame: DataFrame containing patch coordinates for each WSI.
    """
    # Initialize DataFrame to store patch coordinates for all WSIs
    sWindFeatDf = pd.DataFrame()

    # Convert patch shape and stride shape to NumPy arrays for easy manipulation
    patch_shape = np.array(patch_shape)
    stride_shape = np.array(stride_shape)
    
    # Loop through each sample ID
    for sample_id in tqdm(cases, desc="Processing WSIs"):
        # Define WSI path
        wsi_path = f'{WSI_DIR}/{sample_id}.tif'

        WSIReader = OpenSlideWSIReader(wsi_path)
        BASE_RESOLUTION = {"resolution": resolution, "units": "mpp"}
        
        # Get the WSI dimensions at base resolution
        wsi_shape = WSIReader.slide_dimensions(**BASE_RESOLUTION)
        
        # Extract patch coordinates using PatchExtractor
        patch_coords = PatchExtractor.get_coordinates(
            image_shape=wsi_shape,
            patch_input_shape=patch_shape,
            patch_output_shape=None,
            stride_shape=stride_shape,
            input_within_bound=False,
            output_within_bound=True,
        )
        
        # Convert patch coordinates to DataFrame format and add sample ID
        patch_df = pd.DataFrame(patch_coords, columns=['x0', 'y0', 'x1','y1'])
        patch_df['sample_key'] = sample_id
        
        # Append to the main DataFrame
        sWindFeatDf = pd.concat([sWindFeatDf, patch_df], ignore_index=True)
    
    return sWindFeatDf


import geopandas as gpd
from shapely.geometry import box

def get_cells_in_patches(
    cellGeomDf,
    windowDf,
    sample_key_col='sample_key'
    
    ):
    """
    Computes which cells occur in which patch, and returns a DataFrame with patch_id and cell_ids.
    
    Args:
        cellGeomDf (GeoDataFrame): DataFrame containing cell geometries (as polygons) and their coordinates.
        windowDf (DataFrame): DataFrame containing sliding window patch coordinates.
        sample_key (str): Column name in windowDf for unique sample identifier.
        
    Returns:
        DataFrame: DataFrame containing patch_id and cell_ids for each patch.
    """
    samples = list(set(windowDf['sample_key']))
    results = []

    for sample_key in tqdm(samples):
        sampleCellGeomDf = cellGeomDf[cellGeomDf[sample_key_col] == sample_key]
        sampleWindowDf = windowDf[windowDf[sample_key_col] == sample_key]

        # Convert cellGeomDf to a GeoDataFrame if it isn’t already
        if not isinstance(sampleCellGeomDf, gpd.GeoDataFrame):
            sampleCellGeomDf = gpd.GeoDataFrame(sampleCellGeomDf, geometry='geometry')

        # Create a spatial index on cellGeomDf for faster lookups
        sampleCellGeomDf.sindex

        # Process each row in the window DataFrame
        for widx in tqdm(range(sampleWindowDf.shape[0])):
            window = sampleWindowDf.iloc[widx]
            patch_polygon = box(window['x0'], window['y0'], window['x1'], window['y1'])

            # Query the spatial index for cells that may be entirely contained within the patch
            possible_matches_index = list(sampleCellGeomDf.sindex.query(patch_polygon, predicate="contains"))

            # Retrieve the cells that are entirely contained within the patch
            cells_in_patch = sampleCellGeomDf.iloc[possible_matches_index]

            # Create a unique patch_id based on window coordinates and sample_key
            patch_id = f"{sample_key}_{window['x0']}_{window['y0']}_{window['x1']}_{window['y1']}"

            # Add a record for each cell in the patch
            for cell_idx in cells_in_patch.index:
                results.append({
                    'patch_id': patch_id,
                    'cell_id': cell_idx
                })

    # Convert results to DataFrame
    pairDf = pd.DataFrame(results)

    return pairDf


def generate_annotations_patch_level(
    itr_df, 
    feats_cols: List[str], 
    cases: List[str], 
    OUT_ANNOT_DIR: str, 
    scaling: bool = True
):
    """
    Generate annotation files for each sample_id in 'cases', with optional scaling.

    Parameters:
    - itr_df (pd.DataFrame): DataFrame containing indexed data with 'index' column for coordinates and feature columns.
    - feats_cols (List[str]): List of feature columns to include in annotations, along with 'scores'.
    - cases (List[str]): List of sample IDs to filter by in the 'index' column of itr_df.
    - OUT_ANNOT_DIR (str): Directory path where the output database files will be saved.
    - scaling (bool): If True, scales the feature columns between 0 and 1.
    """
    
    os.makedirs(OUT_ANNOT_DIR, exist_ok=True)
    
    # Generate p_coords columns from 'index'
    p_coords = ['x0', 'y0', 'x1', 'y1']
    itr_df[p_coords] = itr_df['index'].str.split('_').str[-4:].apply(pd.Series).astype(int).to_numpy()
    
    # Optional scaling of feature columns
    if scaling:
        itr_df[feats_cols] = MinMaxScaler().fit_transform(itr_df[feats_cols])
        itr_df[feats_cols] = itr_df[feats_cols].clip(upper=1)  # Cap values at 1
    
    # Process each sample_id
    for sample_id in cases:
        annotations = []
        SQ = SQLiteStore()
        
        # Filter DataFrame for the current sample_id
        selDf = itr_df[itr_df['index'].str.contains(sample_id)]
        
        # Create annotations
        for _, row in selDf.iterrows():
            annotations.append(
                Annotation(
                    geometry=Polygon.from_bounds(*row[p_coords].tolist()),
                    properties=row[feats_cols].to_dict()
                )
            )
        
        # Append all annotations and save to file
        SQ.append_many(annotations)
        out_db_file = f'{OUT_ANNOT_DIR}/{sample_id}.db'
        
        # Create index and dump data to file
        SQ.create_index("area", '"area"')
        SQ.dump(out_db_file)

        print(f'Annotation file for sample_id {sample_id} saved to {out_db_file}')


def generate_its_stats(itr_pred_dict, feats_clos):
    rows = []
    for it_id in itr_pred_dict.keys():
        row = {
                'it_id': it_id,
                'score': float(itr_pred_dict[it_id]['scores']),
                'label': itr_pred_dict[it_id]['labels'],
                'fold_idx': itr_pred_dict[it_id]['fold_idx']
            }

        effects = itr_pred_dict[it_id]['effects']
        for feat_name, effect_val in zip(feats_cols, effects):
            row[feat_name] = effect_val
            rows.append(row)
    # Convert to a DataFrame
    df = pd.DataFrame(rows)
    return df