import spatialdata as sd
import os
import pandas as pd
import numpy as np
import geopandas as gpd

save_folder = 'path/to/output'

data_folder = 'path/to/main_csv'

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0)
cluster_df = pd.read_csv('anndata_cc_3_n10.csv', index_col=0)

cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']


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
    cell_types = []
    cluster_type = []
    
    cnt = 0
    
    set_cluster_ids = set(cluster_ids)
    
    in_mask = [name in set_cluster_ids for name in cell_ids_orig]
    
    cell_ids = np.array(cell_ids_orig)[in_mask]
    cell_types = adata.obs["obj.anno_3_w_megs_w_stromal_x"].iloc[in_mask]
    set_cell_ids = set(cell_ids)
    cluster_idx = np.where([name in set_cell_ids for name in cluster_ids])[0]
    cluster_type = cluster_df['cluster_cellcharter_n5'].iloc[cluster_idx]
    
    in_mask = [name in set_cell_ids for name in np.array(cell_centroid.index)]
    cell_centroid_ = cell_centroid[in_mask]
    cell_bound = obj["new_cell_boundaries"][in_mask]
    
    for cc in range(len(cluster_type_list)):
        cluster_to_name = cluster_type_list[cc]
        subsavefolder = '{}/{}'.format(save_folder,cluster_to_name)
        
        if os.path.exists('{}/{}.xlsx'.format(subsavefolder,sample_id)):
            print('{} out of {} : {} --> {} already exists'.format(i+1,len(zarr_),sample_id, cluster_to_name))
            continue
        
        if not os.path.exists(subsavefolder):
            os.makedirs(subsavefolder)
        
        cluster_idx = np.where(np.array(cluster_type) == int(cluster_type_list[cc]))
        cluster_to = cell_bound.iloc[cluster_idx]
        
        cluster_w_cluster = gpd.sjoin_nearest(cell_bound, cluster_to, distance_col="distances", how="left")
        cluster_w_cluster_dropped = cluster_w_cluster.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        cluster_w_cluster_dropped = cluster_w_cluster_dropped.rename(columns={"index_right": "nearest_cluster", "distances": "distance_to_cluster"})
        cluster_w_cluster_dropped = cluster_w_cluster_dropped[["nearest_cluster", "distance_to_cluster"]]
        
        in_mask = [name in set_cell_ids for name in np.array(cluster_w_cluster_dropped['distance_to_cluster'].index)]
        
        distance_cluster = cluster_w_cluster_dropped['distance_to_cluster'][in_mask]
        
        x_data = np.zeros((len(cell_ids),4))
        x_data = x_data.astype(object)
        cluster_type_temp = np.array(cluster_type)
        x_data[:,0] = cluster_type_temp
        x_data[:,1] = cell_centroid_.x
        x_data[:,2] = cell_centroid_.y
        x_data[:,3] = distance_cluster
        header = ['Cluster Type', 'x', 'y', 'Distance to Cluster']
        dff = pd.DataFrame(x_data)
        dff.columns = header
        dff['Distance to Cluster'] = dff['Distance to Cluster'].astype(float)
        dff.set_index('Cluster Type', inplace=True)
        dff.to_excel('{}/{}.xlsx'.format(subsavefolder,sample_id))
