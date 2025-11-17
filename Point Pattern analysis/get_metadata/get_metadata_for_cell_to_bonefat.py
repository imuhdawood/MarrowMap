import spatialdata as sd
import os
import pandas as pd
import numpy as np
import geopandas as gpd

save_folder = 'path/to/output'
data_folder = 'path/to/data'

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0) # Main metadata csv

print('Started PPM')

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    if os.path.exists('{}/{}.xlsx'.format(save_folder,sample_id)):
        print('{} already done'.format(sample_id))
        continue
    
    obj = sd.read_zarr(os.path.join(data_folder,'july_sdata_with_adipocytes', zarr_[i]))
    # df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_20241002_updated_13102024.csv'.format(data_folder), index_col=0)
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
    dist_bone = adata.obs['distance_to_bone']
    fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(sample_id) in x]
    fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]
    
    cell_centroid_x = []
    cell_centroid_y = []
    
    for name, poly in cell_centroid.items():
        if name in cell_ids.values:
            cell_centroid_x.append(poly.x)
            cell_centroid_y.append(poly.y)
    
    cells_w_fat = gpd.sjoin_nearest(cell_bound, fat_gdf, distance_col="distances", how="left")
    fat_dropped = cells_w_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    fat_dropped = fat_dropped.rename(columns={"index_right": "nearest_fat", "distances": "distance_to_fat"})
    fat_dropped = fat_dropped[["nearest_fat", "distance_to_fat"]]
    
    distance_fat = []
    for name, poly in fat_dropped['distance_to_fat'].items():
        if name in cell_ids.values:
            distance_fat.append(poly)
    
    x_data = np.zeros((len(cell_ids),5))
    x_data = x_data.astype(object)
    cell_types_temp = np.array(cell_types)
    x_data[:,0] = cell_types_temp
    x_data[:,1] = cell_centroid_x
    x_data[:,2] = cell_centroid_y
    x_data[:,3] = dist_bone
    x_data[:,4] = distance_fat
    
    header = ['Cell Type', 'x', 'y', 'Distance to Bone', 'Distance to Fat']
    
    dff = pd.DataFrame(x_data)
    dff.columns = header
    dff.set_index('Cell Type', inplace=True)
    dff.to_excel('{}/{}.xlsx'.format(save_folder,sample_id))
    print('{} out of {} : {} has been done'.format(i+1,len(zarr_),sample_id))
