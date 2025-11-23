import spatialdata as sd
import os
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt as dm
import cv2


save_folder = 'path/to.output'
data_folder = 'path/to/main_csv'

cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0)

cluster_df = pd.read_csv('anndata_cc_3_n10.csv', index_col=0)

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    if os.path.exists('{}/{}/{}.csv'.format(save_folder,cluster_type_list[-1],sample_id)):
        print('{} already done'.format(sample_id))
        continue
    
    print('Started {}'.format(sample_id))
    
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
    
    width = int(obj['morphology_focus']['scale0'].dims['x'])
    height = int(obj['morphology_focus']['scale0'].dims['y'])
    
    for ii in range(len(cluster_type_list)):
        
        sub_savefolder = '{}/{}'.format(save_folder,cluster_type_list[ii])
        
        if os.path.exists('{}/{}.csv'.format(sub_savefolder,sample_id)):
            print('{}: {} already exists'.format(sample_id,cluster_type_list[ii]))
            continue
        
        if not os.path.exists(sub_savefolder):
            os.makedirs(sub_savefolder)
        
        img = Image.new('L', (width, height), 0)  # 'L' mode for grayscale image
    
        draw = ImageDraw.Draw(img)
        
        cluster_idx = np.where(np.array(cluster_type)==int(cluster_type_list[ii]))
        
        cell_bound_cluster = cell_bound.iloc[cluster_idx]
        
        for iii in range(len(cell_bound_cluster)):
            try:
                adjusted_polygon = [(x, y) for x, y in cell_bound_cluster['geometry'][iii].exterior.coords]
            except:
                geoms = [g for g in cell_bound_cluster['geometry'][iii].geoms] #List all multipolygon polygons
                maxx = max(geoms, key=lambda x: x.bounds[2])
                adjusted_polygon = [(x, y) for x, y in maxx.exterior.coords]
            draw.polygon(adjusted_polygon, outline=1, fill=1)
        
        binary_mask = np.array(img)
        binary_mask_downsampled = cv2.resize((1-binary_mask),(0,0),fx=1/20,fy=1/20,interpolation=cv2.INTER_NEAREST)
        distance_map = dm(binary_mask_downsampled)
        
        x_coords, y_coords = np.meshgrid(range(1,distance_map.shape[1]+1), range(1,distance_map.shape[0]+1), indexing='xy')
        
        x_flat = x_coords.flatten()
        y_flat = y_coords.flatten()
        
        df_cell = pd.DataFrame({'x': x_flat, 'y': y_flat, 'value': distance_map.flatten()})
        
        df_cell.to_csv('{}/{}.csv'.format(sub_savefolder,sample_id), index=False)
        
        print('{}/{}: {} out of {} has been done'.format(i+1,len(zarr_), ii+1, len(cluster_type_list)))
    
    
    print('{} out of {} has been done'.format(i+1,len(zarr_)))
