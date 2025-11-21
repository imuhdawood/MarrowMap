import spatialdata as sd
import os
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt as dm
import cv2

save_folder = 'path/to/output'
coord_path = 'path/to/ppm'
data_folder = 'path/to/data'

bone_folder = '{}/bone'.format(save_folder)
fat_folder = '{}/fat'.format(save_folder)

if not os.path.exists(save_folder):
    bone_folder = '{}/bone'.format(save_folder)
    fat_folder = '{}/fat'.format(save_folder)
    os.makedirs(bone_folder)
    os.makedirs(fat_folder)

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

for i in range(len(zarr_)):
    if i != 12:
        continue
    sample_id = zarr_[i].split('_with')[0]
    
    print('Started {}'.format(sample_id))
    
    coord_cells = pd.read_excel('{}/{}.xlsx'.format(coord_path,sample_id))    

    obj = sd.read_zarr(os.path.join(data_folder,'july_sdata_with_adipocytes', zarr_[i]))
    polygon_grp = obj['transformed_bone']['geometry']
    
    # pdb.set_trace()
    
    fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(sample_id) in x]
    fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]['geometry']
    
    width = int(obj['morphology_focus']['scale0'].dims['x'])
    height = int(obj['morphology_focus']['scale0'].dims['y'])
    
    del obj
    
    img = Image.new('L', (width, height), 0)  # 'L' mode for grayscale image

    draw = ImageDraw.Draw(img)
    
    for ii in range(len(polygon_grp)):
        try:
            adjusted_polygon = [(x, y) for x, y in polygon_grp[ii].exterior.coords]
        except:
            geoms = [g for g in polygon_grp[ii].geoms] #List all multipolygon polygons
            maxx = max(geoms, key=lambda x: x.bounds[2])
            adjusted_polygon = [(x, y) for x, y in maxx.exterior.coords]
        draw.polygon(adjusted_polygon, outline=1, fill=1)
    
    img_f = Image.new('L', (width, height), 0)  # 'L' mode for grayscale image

    draw_f = ImageDraw.Draw(img_f)
    
    for ii in range(len(fat_gdf)):
        try:
            adjusted_polygon = [(x, y) for x, y in fat_gdf[ii].exterior.coords]
        except:
            geoms = [g for g in fat_gdf[ii].geoms] #List all multipolygon polygons
            maxx = max(geoms, key=lambda x: x.bounds[2])
            adjusted_polygon = [(x, y) for x, y in maxx.exterior.coords]
        draw_f.polygon(adjusted_polygon, outline=1, fill=1)
    
    binary_mask = np.array(img)
    
    binary_mask_f = np.array(img_f)
    
    binary_mask_downsampled = cv2.resize((1-binary_mask),(0,0),fx=1/20,fy=1/20,interpolation=cv2.INTER_NEAREST)
    
    binary_mask_f_downsampled = cv2.resize((1-binary_mask_f),(0,0),fx=1/20,fy=1/20,interpolation=cv2.INTER_NEAREST)
    
    distance_map_f = dm(binary_mask_f_downsampled)
    
    distance_map = dm(binary_mask_downsampled)
    
    pdb.set_trace()
    
    x_coords, y_coords = np.meshgrid(range(1,distance_map.shape[1]+1), range(1,distance_map.shape[0]+1), indexing='xy')
    
    x_flat = x_coords.flatten()
    y_flat = y_coords.flatten()
    
    df = pd.DataFrame({'x': x_flat, 'y': y_flat, 'value': distance_map.flatten()})
    
    df.to_csv('{}/{}.csv'.format(bone_folder,sample_id), index=False)
    
    x_coords_f, y_coords_f = np.meshgrid(range(1,distance_map_f.shape[1]+1), range(1,distance_map_f.shape[0]+1), indexing='xy')
    
    x_flat_f = x_coords_f.flatten()
    y_flat_f = y_coords_f.flatten()
    
    df_f = pd.DataFrame({'x': x_flat_f, 'y': y_flat_f, 'value': distance_map_f.flatten()})
    
    df_f.to_csv('{}/{}.csv'.format(fat_folder,sample_id), index=False)
    
    print('{} out of {} has been done'.format(i+1,len(zarr_)))
