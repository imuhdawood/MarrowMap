import os
import spatialdata as sd
import numpy as np
import pandas as pd
from shapely.geometry import Point

score_path = 'path/to/cif_score'

distance_path = 'path/to/distance_map/bone'

transform_path = 'path/to/landmarks_matrix'

data_folder = 'path/to/main_data'

stride_ratio = 4

sub_tile_size = (512/(stride_ratio))*(1/20)

min_dist = (((sub_tile_size)**2)*2)**(1/2)

save_path = 'path/to/output'

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

if not os.path.exists(save_path):
    os.makedirs(save_path)


dist_folder_list = []

for file in os.listdir(distance_path):
    dist_folder_list.append(file)

df_header = ['Distance', 'CIF_Score', 'ITS']

for i in range(len(dist_folder_list)):
    upper_id = dist_folder_list[i].split('_')[0]
    lower_id = dist_folder_list[i].split('.')[0]
    
    if os.path.exists('{}/{}.csv'.format(save_path,lower_id)):
        print('{} already done'.format(lower_id))
        continue
    
    print('{} is started'.format(lower_id))
    
    with open('{}/{}/cif_scores.txt'.format(score_path,upper_id)) as p:
        strings_p = p.readlines()
        score_x  = [int(string.split(' ')[0].split(':')[-1]) for string in strings_p[:-1]]
        score_y  = [int(string.split(' ')[1].split(':')[-1]) for string in strings_p[:-1]]
        tile_score = [float(string.split(' ')[2].split(':')[-1]) for string in strings_p[:-1]]
        p.close()
    
    obj = sd.read_zarr(os.path.join(data_folder,'july_sdata_with_adipocytes', '{}_with_it_cellularity.zarr'.format(lower_id)))
    
    fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(lower_id) in x]
    fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]
    fat_poly_grp = fat_gdf['geometry']
    
    del fat_cells
    del fat_gdf
    
    its_poly_grp = obj['intertrabecular_regions']['geometry']
    
    del obj
    
    trnsmatrix = np.array(pd.read_csv('{}/{}_matrix.csv'.format(transform_path,lower_id),header=None))
    
    inv_trns = np.linalg.inv(trnsmatrix)
    
    
    dist_xl = pd.read_csv('{}/{}'.format(distance_path,dist_folder_list[i]))
    max_dist = max(dist_xl['value'])
    dist_cnt = int(np.floor(max_dist/min_dist)+1)
    
    aa = np.zeros((max(dist_xl['y']),max(dist_xl['x'])))
    for aaa in range(len(dist_xl['x'])):
        aa[dist_xl['y'][aaa]-1,dist_xl['x'][aaa]-1] = dist_xl['value'][aaa]
    
    anch = 0
    anch_sub = 0
    
    df_temp = np.zeros((1,3))
    anch_df = 0
    
    min_dist_cnt = 0
    
    for ii in range(len(score_x)):
        score_coord_mat = [int(score_x[ii]+512/(stride_ratio*2)),int(score_y[ii]+512/(stride_ratio*2)),1]
        
        trns_coord_x, trns_coord_y, _ = np.dot(trnsmatrix, score_coord_mat)
        
        try:
            dist_tile = np.array(dist_xl['value'])[np.where((np.array(dist_xl['x'])==int(trns_coord_x/20))&(np.array(dist_xl['y'])==int(trns_coord_y/20)))[0][0]]
        except:
            print('{}/{} {}: {} out of {} has been done'.format(i+1,len(dist_folder_list), lower_id, ii+1, len(score_x)))
            continue
        
        if dist_tile <= ((sub_tile_size/2)**2*2)**(1/2):
            print('{}/{} {}: {} out of {} has been done'.format(i+1,len(dist_folder_list), lower_id, ii+1, len(score_x)))
            continue
        
        which_its=[]
        in_adipo = 0
        
        point = Point(trns_coord_x, trns_coord_y)
        
        for iv in range(len(its_poly_grp)):
            if its_poly_grp[iv].contains(point):
                # pdb.set_trace()
                which_its = int(its_poly_grp.index[iv].split('_R')[-1])
                continue
        
        if not which_its:
            print('{}/{} {}: {} out of {} has been done'.format(i+1,len(dist_folder_list), lower_id, ii+1, len(score_x)))
            continue
        
        for iv in range(len(fat_poly_grp)):
            if fat_poly_grp[iv].contains(point):
                # pdb.set_trace()
                in_adipo = 1
                continue
        
        if in_adipo == 1:
            print('{}/{} {}: {} out of {} has been done'.format(i+1,len(dist_folder_list), lower_id, ii+1, len(score_x)))
            continue
        
        
        if anch_df == 0:        
            df_temp[0,0] = dist_tile
            # pdb.set_trace()
            df_temp[0,1] = tile_score[ii]
            df_temp[0,2] = which_its
            anch_df = 1
        else:
            df_tempp = np.zeros((1,3))
            df_tempp[0,0] = dist_tile
            df_tempp[0,1] = tile_score[ii]
            df_tempp[0,2] = which_its
            df_temp = np.concatenate((df_temp,df_tempp))
        
        print('{}/{} {}: {} out of {} has been done'.format(i+1,len(dist_folder_list), lower_id, ii+1, len(score_x)))
        
    
    df = pd.DataFrame(df_temp)
    df.columns = df_header
    df.to_csv('{}/{}.csv'.format(save_path,lower_id), index=False)
    print('{} out of {} has been done'.format(i+1,len(dist_folder_list)))
