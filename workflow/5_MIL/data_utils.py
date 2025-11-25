import numpy as np
import pandas as pd
from tqdm import tqdm
import os
from shapely import wkt

def create_bags_from_dataframe(df, feats_cols, labels_col='Normal_vs_MPN', bag_id='sample_key'):
   
    # Group the data by bag_id (sample_key)
    grouped = df.groupby(bag_id)

    # Initialize a dictionary to store bag-related information
    bag_data = {
        'samples_id': [], 
        'bags': [], 
        'labels': [], 
        'samples_its_id': [], 
        'samples_its_labels': []
    }
    
    # Iterate over each group to extract features and labels
    for name, group in grouped:
        # Add features (bag), label, and other metadata to the dictionary directly
        bag_data['samples_id'].append(name)
        bag_data['bags'].append(group[feats_cols].to_numpy())
        bag_data['labels'].append(group[labels_col].iloc[0])
        bag_data['samples_its_id'].append(group.index.tolist())
        bag_data['samples_its_labels'].append(group[labels_col].tolist())

    # Convert the list values to numpy arrays
    bag_data = {key: np.array(value, dtype='object') for key, value in bag_data.items()}

    return bag_data

def process_intertrabecular_data(
    ITS_POLY_FILE, 
    INPUT_FILE, 
    OUT_FILE, 
    rkey, 
    skey, 
    dkey, 
    ENVIRON, 
    CELL_TYPE, 
    poly_to_area,
    sample_itr_weight_map=None,
    environ_selected_feats = None
    ):
    """
    Processes intertrabecular data and computes IT-level properties.
    
    Args:
        ITS_POLY_FILE (str): Path to the polygon geometry file.
        INPUT_FILE (str): Path to the input file containing sample data.
        OUT_FILE (str): Path to save the output file. If a file with this name already exists, it will be loaded instead of recalculating.
        rkey (str): Column name for region key in the data.
        skey (str): Column name for sample key in the data.
        dkey (str): Column name for additional data key.
        ENVIRON (str): Column name for environment property.
        CELL_TYPE (str): Column name for cell type property.
        poly_to_area (function): Function to compute the area of a polygon.
        sample_itr_weight_map (dict, optional): Dictionary to store sample weights. Pass an empty dict to populate weights.

    Returns:
        itrDf (DataFrame): Processed DataFrame with computed IT-level properties.
    """
    # Check if the output file already exists, load it if it does
    if os.path.isfile(OUT_FILE):
        print("Loading existing output file.")
        return pd.read_csv(OUT_FILE)

    try:
        # Load geometry and data files
        itrGeom = pd.read_csv(ITS_POLY_FILE, index_col='Unnamed: 0')
        consolDf = pd.read_csv(INPUT_FILE, index_col='Unnamed: 0')
        consolDf[skey] = consolDf['cell_id'].str.split('_').str[-2:].str.join('_')
        # Generate sample keys and filter data
        consolDf = consolDf[~consolDf[rkey].isin(['non_intertrabecular'])]
        consolDf.index = consolDf[skey].astype(str) + '_R' + consolDf[rkey].astype(float).astype(int).astype(str)
        consolDf[consolDf[CELL_TYPE].reset_index().groupby('index')[CELL_TYPE].count()>=100]
        # Selecting Interesting Miceoenviron Types

        if environ_selected_feats:
            consolDf = consolDf[consolDf[ENVIRON].isin(environ_selected_feats)]
        
        # Join geometry data and remove duplicates
        itrDf = consolDf[[rkey, skey, dkey]].drop_duplicates().join(itrGeom).dropna()

        # Identify unique sample keys
        samples_keys = itrDf.index.str.split('_').str[:2].str.join('_').unique()
       
        # Calculate weights if sample_itr_weight_map is provided
        if sample_itr_weight_map is not None:
            sample_itr_weight = pd.DataFrame(index=itrGeom.index)
            sample_itr_weight['area'] = itrGeom['geometry'].apply(poly_to_area)
            
            # Normalize weights by sample
            sample_itr_weight_map.clear()  # Clear any existing data
            for skey in samples_keys:
                subDf = sample_itr_weight[sample_itr_weight.index.str.contains(skey)]
                sample_itr_weight_map.update((subDf / subDf.sum())['area'].to_dict())
        
        # Compute IT-level properties
        for key in tqdm(itrDf.index):
            ITDf = consolDf.loc[[key], :].dropna(subset=[ENVIRON, CELL_TYPE])
            
            # Get weighted or unweighted counts based on sample_itr_weight_map availability
            if sample_itr_weight_map is not None:
                weight = sample_itr_weight_map[key]
                type_counts = ITDf[ENVIRON].value_counts(normalize=True) * weight
                cell_counts = ITDf[CELL_TYPE].value_counts(normalize=True) * weight
            else:
                type_counts = ITDf[ENVIRON].value_counts(normalize=True)
                cell_counts = ITDf[CELL_TYPE].value_counts(normalize=True)
            
            # Assign values to the main DataFrame
            itrDf.loc[key, type_counts.index] = type_counts.values
            itrDf.loc[key, cell_counts.index] = cell_counts.values
        
        # Fill NaN values with 0 and save to the output file
        itrDf.fillna(0, inplace=True)
        itrDf.to_csv(OUT_FILE, index=None)
        print("ITS level stats saved successfully.")
        
        return itrDf
    
    except Exception as e:
        print(f"An error occurred while reading file: {e}")
        return None
    
def compute_patch_stats(
    cases,
    WSI_DIR,
    CELL_GEOM_FILE,
    INPUT_FILE,
    ENVIRON,
    CELL_TYPE,
    PATCH_OUT_FILE,
    extract_patches_from_wsi,
    get_cells_in_patches, 
    patch_size = (1024, 1024),
    stride = (512, 512)
):
    """
    Compute patch-level statistics for cells in Whole Slide Images (WSI) based on cell geometries and spatial windows.

    Args:
        cases (list): List of WSI sample IDs to process.
        WSI_DIR (str): Directory containing WSI files.
        CELL_GEOM_FILE (str): Path to the CSV file containing cell geometries.
        INPUT_FILE (str): Path to the CSV file containing additional cell data.
        ENVIRON (str): Column name for the cell environment type.
        CELL_TYPE (str): Column name for the cell type.
        PATCH_OUT_FILE (str): Path to save the output CSV with patch-level statistics.
        extract_patches_from_wsi (function): Function to extract spatial windows from WSI.
        get_cells_in_patches (function): Function to map cells to patches based on geometry.

    Returns:
        None: Saves computed patch-level statistics to the specified output file.
    """
    # Check if output file exists and load if available
    if os.path.isfile(PATCH_OUT_FILE):
        print('Loading existing patch statistics...')
        patchStatsDf = pd.read_csv(PATCH_OUT_FILE, index_col = 'Unnamed: 0')
        return patchStatsDf
    else:
        # Extract patches/windows from WSI and load cell geometries
        sWindowCoordDf = extract_patches_from_wsi(
            cases,
            WSI_DIR,
            patch_shape = patch_size, 
            stride_shape = stride
        )
        cellGeom = pd.read_csv(CELL_GEOM_FILE, index_col='Unnamed: 0')
        
        columns = ['x', 'y', 'geometry', 'sample_key']
        cellGeom = cellGeom[columns].dropna()
        cellGeom['geometry'] = cellGeom['geometry'].apply(wkt.loads)
        consolDf = pd.read_csv(INPUT_FILE).set_index('cell_id').join(cellGeom).dropna(subset=[ENVIRON, CELL_TYPE])
        
        # Map cells to patches and filter by unique cell count threshold
        windowCellsDf = get_cells_in_patches(consolDf[columns], sWindowCoordDf)
        windowCellsDf = windowCellsDf[windowCellsDf.groupby('patch_id')['cell_id'].transform('nunique') > 50]
        windowCellsDf.set_index('patch_id', inplace=True)

        # Initialize DataFrame for patch-level statistics
        patchStatsDf = pd.DataFrame()

        # Process patches to compute environment and cell type proportions
        for patch_id, patch_data in tqdm(windowCellsDf.groupby(level=0), desc="Processing patches"):
            patch_cells = patch_data['cell_id'].tolist()
            patchDf = consolDf.loc[patch_cells]

            # Calculate normalized counts for ENVIRON and CELL_TYPE columns
            type_counts = patchDf[ENVIRON].value_counts(normalize=True)
            cell_counts = patchDf[CELL_TYPE].value_counts(normalize=True)

            # Store statistics in patchStatsDf
            patchStatsDf.loc[patch_id, type_counts.index] = type_counts.values
            patchStatsDf.loc[patch_id, cell_counts.index] = cell_counts.values

        # Fill NaN values with 0, cases where a certain feature is not present in a patch
        patchStatsDf.fillna(0, inplace=True)
        patchStatsDf.to_csv(PATCH_OUT_FILE)
        print(f"Patch-level statistics saved to {PATCH_OUT_FILE}")
        return patchStatsDf