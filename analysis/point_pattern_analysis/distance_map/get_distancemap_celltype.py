import spatialdata as sd
import os
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt as dm
import cv2

save_folder = 'path/to/output'
coord_path = 'path/to/ppm' # From get_metadata
data_folder = 'path/to/data'

"""
cell_type:
    
The default is 'All' which generates distance maps of all cell types
If you just want to generate distanc map for specific cell type, 
copy one of the cell type from cell_types_list
"""
cell_type = 'All'

cell_types_list = ['Adipo-MSC','B_cell', 'DC', 'Endothelial', 'Erythroid', 'GMP','Granulocyte/mast', 'HSPC',
                   'Macrophage', 'Megakaryocyte','MNP', 'Monocyte', 'Myeloid', 'Osteo-MSC', 'Plasma_cell', 
                   'SMC', 'Stromal', 'T_cell']

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0)

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    if os.path.exists('{}/{}/{}.csv'.format(save_folder,cell_types_list[-1],sample_id)):
        print('{} already done'.format(sample_id))
        continue
    
    print('Started {}'.format(sample_id))
    
    coord_cells = pd.read_excel('{}/{}.xlsx'.format(coord_path,sample_id))    

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
    
    if cell_type == 'All':
        for ii in range(len(cell_types_list)):
            if '/' not in cell_types_list[ii]:
                sub_savefolder = '{}/{}'.format(save_folder,cell_types_list[ii])
            else:
                sub_savefolder = '{}/{}_{}'.format(save_folder,cell_types_list[ii].split('/')[0],
                                                   cell_types_list[ii].split('/')[1])
            
            if os.path.exists('{}/{}.csv'.format(sub_savefolder,sample_id)):
                print('{}: {} already exists'.format(sample_id,cell_types_list[ii]))
                continue
            
            if not os.path.exists(sub_savefolder):
                os.makedirs(sub_savefolder)
            
            cell_boundaries = []
            cnt = 0
            for name, poly in cell_bound['geometry'].items():
                if name in cell_ids.values:
                    if obj["table"].obs["obj.anno_3_w_megs_w_stromal_x"][cnt] == cell_types_list[ii]:
                        cell_boundaries.append(poly)
                        cnt+=1
                    else:
                        cnt+=1
                else:
                    cnt+=1
        
            width = int(obj['morphology_focus']['scale0'].dims['x'])
            height = int(obj['morphology_focus']['scale0'].dims['y'])
            
            img = Image.new('L', (width, height), 0)  # 'L' mode for grayscale image
        
            draw = ImageDraw.Draw(img)
            
            for iii in range(len(cell_boundaries)):
                try:
                    adjusted_polygon = [(x, y) for x, y in cell_boundaries[iii].exterior.coords]
                except:
                    geoms = [g for g in cell_boundaries[iii].geoms] #List all multipolygon polygons
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
            
            print('{}/{}: {} out of {} has been done'.format(i+1,len(zarr_), ii+1, len(cell_types_list)))
            
    else:
        if '/' not in cell_type:
            sub_savefolder = '{}/{}'.format(save_folder,cell_type)
        else:
            sub_savefolder = '{}/{}_{}'.format(save_folder,cell_type.split('/')[0],
                                               cell_type.split('/')[1])
        
        if os.path.exists('{}/{}.csv'.format(sub_savefolder,sample_id)):
            print('{}: {} already exists'.format(sample_id, cell_type))
            continue
        
        if not os.path.exists(sub_savefolder):
            os.makedirs(sub_savefolder)
        
        cell_boundaries = []
        cnt = 0
        for name, poly in cell_bound['geometry'].items():
            if name in cell_ids.values:
                if obj["table"].obs["obj.anno_3_w_megs_w_stromal_x"][cnt] == cell_type:
                    cell_boundaries.append(poly)
                    cnt+=1
                else:
                    cnt+=1
            else:
                cnt+=1
    
        width = int(obj['morphology_focus']['scale0'].dims['x'])
        height = int(obj['morphology_focus']['scale0'].dims['y'])
        
        img = Image.new('L', (width, height), 0)  # 'L' mode for grayscale image
    
        draw = ImageDraw.Draw(img)
        
        for iii in range(len(cell_boundaries)):
            try:
                adjusted_polygon = [(x, y) for x, y in cell_boundaries[iii].exterior.coords]
            except:
                geoms = [g for g in cell_boundaries[iii].geoms] #List all multipolygon polygons
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
        
        print('{}/{}: {} : has been done'.format(i+1,len(zarr_), cell_type))
        
    print('{} out of {} has been done'.format(i+1,len(zarr_)))
