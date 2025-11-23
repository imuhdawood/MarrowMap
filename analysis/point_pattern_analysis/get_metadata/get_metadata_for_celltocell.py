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

cell_list = ['Adipo-MSC','B_cell', 'DC', 'Endothelial', 'Erythroid', 'GMP','Granulocyte/mast', 'HSPC',
             'Macrophage', 'Megakaryocyte', 'Monocyte', 'Myeloid', 'Osteo-MSC', 'Plasma_cell', 
             'SMC', 'Stromal', 'T_cell']

print('Started PPM')

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    if os.path.exists('{}/T_cell/{}.xlsx'.format(save_folder,sample_id)):
        print('{} already done'.format(sample_id))
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
    
    cell_bound_aft = []
    cnt_i = 0
   
    cnt_p = 0
    
    cell_bound_aft = gpd.GeoDataFrame(index = range(1),columns=['geometry'], geometry='geometry')
    
    for name, poly in cell_bound['geometry'].items():
        if name in cell_ids.values:
            if cnt_p == 0:
                cell_bound_aft.iloc[[0]] =cell_bound.iloc[[cnt_i]]
                cell_bound_aft.index = [cell_bound.iloc[[cnt_i]].index[0]]
                cnt_p += 1
            else:
                cell_bound_aft = pd.concat([cell_bound_aft,cell_bound.iloc[[cnt_i]]])
            cnt_i += 1
        else:
            cnt_i += 1
    
    for cc in range(len(cell_list)):
        cnt_i = 0
       
        cnt_p = 0
        cell_to_name = cell_list[cc]
        subsavefolder = '{}/{}'.format(save_folder,cell_to_name)
        
        if os.path.exists('{}/{}.xlsx'.format(subsavefolder,sample_id)):
            print('{} out of {} : {} --> {} already exists'.format(i+1,len(zarr_),sample_id, cell_to_name))
            continue
        
        if not os.path.exists(subsavefolder):
            os.makedirs(subsavefolder)
        
        cell_to = gpd.GeoDataFrame(index = range(1),columns=['geometry'], geometry='geometry')
        for ccc in range(len(cell_types)):
            if cell_types[ccc] == cell_to_name:
                if cnt_p == 0:
                    cell_to.iloc[[0]] =cell_bound_aft.iloc[[cnt_i]]
                    cell_to.index = [cell_bound_aft.iloc[[cnt_i]].index[0]]
                    cnt_p += 1
                else:
                    cell_to = pd.concat([cell_to,cell_bound_aft.iloc[[cnt_i]]])
                cnt_i += 1
            else:
                cnt_i += 1
        
        cells_w_cell = gpd.sjoin_nearest(cell_bound_aft, cell_to, distance_col="distances", how="left")
        cell_dropped = cells_w_cell.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        cell_dropped = cell_dropped.rename(columns={"index_right": "nearest_cell", "distances": "distance_to_cell"})
        cell_dropped = cell_dropped[["nearest_cell", "distance_to_cell"]]
        distance_cell = []
        for name, poly in cell_dropped['distance_to_cell'].items():
            if name in cell_ids.values:
                distance_cell.append(poly)
        
        x_data = np.zeros((len(cell_ids),2))
        x_data = x_data.astype(object)
        cell_types_temp = np.array(cell_types)
        x_data[:,0] = cell_types_temp
        x_data[:,1] = distance_cell
        header = ['Cell Type', 'Distance to Cell']
        dff = pd.DataFrame(x_data)
        dff.columns = header
        dff['Distance to Cell'] = dff['Distance to Cell'].astype(float)
        dff.set_index('Cell Type', inplace=True)
        dff.to_excel('{}/{}.xlsx'.format(subsavefolder,sample_id))
        print('{} out of {} : {} --> {} has been done'.format(i+1,len(zarr_),sample_id, cell_to_name))
