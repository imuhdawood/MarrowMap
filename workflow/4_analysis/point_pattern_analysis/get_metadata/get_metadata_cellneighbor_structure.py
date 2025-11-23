import spatialdata as sd
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely

art_save_folder = 'path/to/arteriole_output'
sinu_save_folder = 'path/to/sinusoid_output'
fat_save_folder = 'path/to/fat_output'
bone_save_folder = 'path/to/bone_output'

data_folder = 'path/to/input_data'
art_folder = 'path/to/arteriole'
sinu_folder = 'path/to/sinosoids'

if not os.path.exists(art_save_folder):
    os.makedirs(art_save_folder)

if not os.path.exists(sinu_save_folder):
    os.makedirs(sinu_save_folder)

if not os.path.exists(fat_save_folder):
    os.makedirs(fat_save_folder)

if not os.path.exists(bone_save_folder):
    os.makedirs(bone_save_folder)

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0)
cluster_df = pd.read_csv('anndata_cc_3_n10.csv', index_col=0)

print('Started PPM')

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    
    obj = sd.read_zarr(os.path.join(data_folder,'july_sdata_with_adipocytes', zarr_[i]))
    adata = obj["table"]
    adata.obs = adata.obs.merge(df, how="left", left_index=True, right_index=True)
    
    df_temp = df[["extra_new_cell_status", "obj.anno_3_w_megs_w_stromal"]]
     
    adata.obs = adata.obs.merge(df_temp, how="left", left_index=True, right_index=True)
    
    adata = adata[adata.obs["annotation_x"] != "remove"]
    adata = adata[adata.obs["annotation_y"] != "remove"]
    adata = adata[adata.obs["extra_new_cell_status_x"] != "artefact"]
    adata = adata[adata.obs["extra_new_cell_status_y"] != "artefact"]
    adata = adata[adata.obs["obj.anno_3_w_megs_w_stromal_x"] != "remove"]
    adata = adata[adata.obs["obj.anno_3_w_megs_w_stromal_y"] != "remove"]
    adata = adata[adata.obs["obj.anno_3_w_megs_w_stromal_x"].isna() != True]
    adata = adata[adata.obs["obj.anno_3_w_megs_w_stromal_y"].isna() != True]
    
    cell_centroid = obj["new_cell_boundaries"]['geometry'].centroid
    cell_ids_orig = adata.obs["cell_id_x"]
    
    cluster_ids = np.array(cluster_df.index)
    
    cell_ids = []
    cell_region_ids = []
    cell_types = []
    dist_bone = []
    cluster_type = []
    
    cnt = 0
    
    set_cluster_ids = set(cluster_ids)
    
    in_mask = [name in set_cluster_ids for name in cell_ids_orig]
    
    cell_ids = np.array(cell_ids_orig)[in_mask]
    cell_types = adata.obs["obj.anno_3_w_megs_w_stromal_x"].iloc[in_mask]
    dist_bone = adata.obs['distance_to_bone'].iloc[in_mask]
    
    set_cell_ids = set(cell_ids)
    cluster_idx = np.where([name in set_cell_ids for name in cluster_ids])[0]
    cluster_type = cluster_df['cluster_cellcharter_n5'].iloc[cluster_idx]
    cell_bound = obj["new_cell_boundaries"]
    fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(sample_id) in x]
    fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]
    
    in_mask = [name in set_cell_ids for name in np.array(cell_centroid.index)]
    cell_centroid_ = cell_centroid[in_mask]
    
    x_data = np.zeros((len(cell_ids),6))
    x_data = x_data.astype(object)
    cell_types_temp = np.array(cell_types)
    x_data[:,0] = cell_ids
    x_data[:,1] = cell_types_temp
    x_data[:,2] = np.array(cluster_type)
    x_data[:,3] = np.array(cell_centroid_.x)
    x_data[:,4] = np.array(cell_centroid_.y)
    x_data[:,5] = np.array(dist_bone)
    
    header = ['Cell ID','Cell Type', 'Cluster Type', 'x', 'y', 'Distance to Bone']
    
    dff = pd.DataFrame(x_data)
    dff.columns = header
    dff.set_index('Cell ID', inplace=True)
    dff.to_excel('{}/{}.xlsx'.format(bone_save_folder,sample_id))
    
    cells_w_fat = gpd.sjoin_nearest(cell_bound, fat_gdf, distance_col="distances", how="left")
    fat_dropped = cells_w_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    fat_dropped = fat_dropped.rename(columns={"index_right": "nearest_fat", "distances": "distance_to_fat"})
    fat_dropped = fat_dropped[["nearest_fat", "distance_to_fat"]]
    
    in_mask = [name in set_cell_ids for name in np.array(fat_dropped['distance_to_fat'].index)]
    
    distance_fat = fat_dropped['distance_to_fat'][in_mask]
    
    x_data = np.zeros((len(cell_ids),6))
    x_data = x_data.astype(object)
    cell_types_temp = np.array(cell_types)
    x_data[:,0] = cell_ids
    x_data[:,1] = cell_types_temp
    x_data[:,2] = np.array(cluster_type)
    x_data[:,3] = np.array(cell_centroid_.x)
    x_data[:,4] = np.array(cell_centroid_.y)
    x_data[:,5] = np.array(distance_fat)
    
    header = ['Cell ID','Cell Type', 'Cluster Type', 'x', 'y', 'Distance to Fat']
    
    dff = pd.DataFrame(x_data)
    dff.columns = header
    dff.set_index('Cell ID', inplace=True)
    dff.to_excel('{}/{}.xlsx'.format(fat_save_folder,sample_id))
    
    if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
    
        art_gdf_temp = pd.read_csv('{}/{}_A.csv'.format(art_folder,sample_id))
    
        art_gdf = gpd.GeoDataFrame(art_gdf_temp, geometry=shapely.wkt.loads(art_gdf_temp['geometry']))
        
        cells_w_art = gpd.sjoin_nearest(cell_bound, art_gdf, distance_col="distances", how="left")
        art_dropped = cells_w_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        art_dropped = art_dropped.rename(columns={"index_right": "nearest_arteriole", "distances": "distance_to_arteriole"})
        art_dropped = art_dropped[["nearest_arteriole", "distance_to_arteriole"]]
        
        in_mask = [name in set_cell_ids for name in np.array(art_dropped['distance_to_arteriole'].index)]
        
        distance_art = art_dropped['distance_to_arteriole'][in_mask]
        
        x_data = np.zeros((len(cell_ids),6))
        x_data = x_data.astype(object)
        cell_types_temp = np.array(cell_types)
        x_data[:,0] = cell_ids
        x_data[:,1] = cell_types_temp
        x_data[:,2] = np.array(cluster_type)
        x_data[:,3] = np.array(cell_centroid_.x)
        x_data[:,4] = np.array(cell_centroid_.y)
        x_data[:,5] = np.array(distance_art)
        
        header = ['Cell ID','Cell Type', 'Cluster Type', 'x', 'y', 'Distance to Arteriole']
        
        dff = pd.DataFrame(x_data)
        dff.columns = header
        dff.set_index('Cell ID', inplace=True)
        dff.to_excel('{}/{}.xlsx'.format(art_save_folder,sample_id))
        print('{} out of {} : {}_A has been done'.format(i+1,len(zarr_),sample_id))
    
    if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
        sinu_gdf_temp = pd.read_csv('{}/{}_S.csv'.format(sinu_folder,sample_id))
    
        sinu_gdf = gpd.GeoDataFrame(sinu_gdf_temp, geometry=shapely.wkt.loads(sinu_gdf_temp['geometry']))
        
        cells_w_sinu = gpd.sjoin_nearest(cell_bound, sinu_gdf, distance_col="distances", how="left")
        sinu_dropped = cells_w_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        sinu_dropped = sinu_dropped.rename(columns={"index_right": "nearest_sinusoid", "distances": "distance_to_sinusoid"})
        sinu_dropped = sinu_dropped[["nearest_sinusoid", "distance_to_sinusoid"]]
        
        in_mask = [name in set_cell_ids for name in np.array(sinu_dropped['distance_to_sinusoid'].index)]
        
        distance_sinu = sinu_dropped['distance_to_sinusoid'][in_mask]
        
        x_data = np.zeros((len(cell_ids),6))
        x_data = x_data.astype(object)
        cell_types_temp = np.array(cell_types)
        x_data[:,0] = cell_ids
        x_data[:,1] = cell_types_temp
        x_data[:,2] = np.array(cluster_type)
        x_data[:,3] = np.array(cell_centroid_.x)
        x_data[:,4] = np.array(cell_centroid_.y)
        x_data[:,5] = np.array(distance_sinu)
        
        header = ['Cell ID','Cell Type', 'Cluster Type', 'x', 'y', 'Distance to Sinusoid']
        
        dff = pd.DataFrame(x_data)
        dff.columns = header
        dff.set_index('Cell ID', inplace=True)
        dff.to_excel('{}/{}.xlsx'.format(sinu_save_folder,sample_id))
        print('{} out of {} : {}_S has been done'.format(i+1,len(zarr_),sample_id))
