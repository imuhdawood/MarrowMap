import spatialdata as sd
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely

art_save_folder = 'path/to/arteriole_output'
sinu_save_folder = 'path/to/sinusoid_output'
data_folder = 'path/to/input_data'
art_folder = 'path/to/arteriole_input' # Transformed arteriole
sinu_folder = 'path/to/sinusoid_input' # Transformed sinusoid

if not os.path.exists(art_save_folder):
    os.makedirs(art_save_folder)

if not os.path.exists(sinu_save_folder):
    os.makedirs(sinu_save_folder)

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0)

print('Started PPM')

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    
    if not os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
        if not os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
            continue
    
    if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
        if os.path.exists('{}/{}.xlsx'.format(art_save_folder,sample_id)):
            print('{}_A already done'.format(sample_id))
            if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
                if os.path.exists('{}/{}.xlsx'.format(sinu_save_folder,sample_id)):
                    print('{}_S already done'.format(sample_id))
                    continue
            else:
                continue
    
    elif os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
        if os.path.exists('{}/{}.xlsx'.format(sinu_save_folder,sample_id)):
            print('{}_S already done'.format(sample_id))
            continue
    
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
    cell_ids = adata.obs["cell_id_x"]
    cell_types = adata.obs["obj.anno_3_w_megs_w_stromal_x"]
    cell_bound = obj["new_cell_boundaries"]
    cell_centroid_x = []
    cell_centroid_y = []
    
    for name, poly in cell_centroid.items():
        if name in cell_ids.values:
            cell_centroid_x.append(poly.x)
            cell_centroid_y.append(poly.y)
    
    if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
    
        art_gdf_temp = pd.read_csv('{}/{}_A.csv'.format(art_folder,sample_id))
    
        art_gdf = gpd.GeoDataFrame(art_gdf_temp, geometry=shapely.wkt.loads(art_gdf_temp['geometry']))
        
        cells_w_art = gpd.sjoin_nearest(cell_bound, art_gdf, distance_col="distances", how="left")
        art_dropped = cells_w_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        art_dropped = art_dropped.rename(columns={"index_right": "nearest_arteriole", "distances": "distance_to_arteriole"})
        art_dropped = art_dropped[["nearest_arteriole", "distance_to_arteriole"]]
        
        distance_art = []
        for name, poly in art_dropped['distance_to_arteriole'].items():
            if name in cell_ids.values:
                distance_art.append(poly)
        
        x_data = np.zeros((len(cell_ids),4))
        x_data = x_data.astype(object)
        cell_types_temp = np.array(cell_types)
        x_data[:,0] = cell_types_temp
        x_data[:,1] = cell_centroid_x
        x_data[:,2] = cell_centroid_y
        x_data[:,3] = distance_art
        
        header = ['Cell Type', 'x', 'y', 'Distance to Arteriole']
        
        dff = pd.DataFrame(x_data)
        dff.columns = header
        dff.set_index('Cell Type', inplace=True)
        dff.to_excel('{}/{}.xlsx'.format(art_save_folder,sample_id))
        print('{} out of {} : {}_A has been done'.format(i+1,len(zarr_),sample_id))
    
    if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
        sinu_gdf_temp = pd.read_csv('{}/{}_S.csv'.format(sinu_folder,sample_id))
    
        sinu_gdf = gpd.GeoDataFrame(sinu_gdf_temp, geometry=shapely.wkt.loads(sinu_gdf_temp['geometry']))
        
        cells_w_sinu = gpd.sjoin_nearest(cell_bound, sinu_gdf, distance_col="distances", how="left")
        sinu_dropped = cells_w_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        sinu_dropped = sinu_dropped.rename(columns={"index_right": "nearest_sinusoid", "distances": "distance_to_sinusoid"})
        sinu_dropped = sinu_dropped[["nearest_sinusoid", "distance_to_sinusoid"]]
        
        distance_sinu = []
        for name, poly in sinu_dropped['distance_to_sinusoid'].items():
            if name in cell_ids.values:
                distance_sinu.append(poly)
        
        x_data = np.zeros((len(cell_ids),4))
        x_data = x_data.astype(object)
        cell_types_temp = np.array(cell_types)
        x_data[:,0] = cell_types_temp
        x_data[:,1] = cell_centroid_x
        x_data[:,2] = cell_centroid_y
        x_data[:,3] = distance_sinu
        
        header = ['Cell Type', 'x', 'y', 'Distance to Sinusoid']
        
        dff = pd.DataFrame(x_data)
        dff.columns = header
        dff.set_index('Cell Type', inplace=True)
        dff.to_excel('{}/{}.xlsx'.format(sinu_save_folder,sample_id))
        print('{} out of {} : {}_S has been done'.format(i+1,len(zarr_),sample_id))
