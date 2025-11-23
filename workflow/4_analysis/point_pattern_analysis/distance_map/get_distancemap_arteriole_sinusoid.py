import os
import spatialdata as sd
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt as dm
import shapely.wkt as wkt
import cv2

structure = ['sinusoid','arteriole']

data_folder = 'path/to/main_csv'
transformed_structure_folder = 'path/to/transformed_structure'
art_save_folder = 'path/to/output_arteriole'
sinu_save_folder = 'path/to/output_sinusoid'

if not os.path.exists(art_save_folder):
    os.makedirs(art_save_folder)
    
if not os.path.exists(sinu_save_folder):
    os.makedirs(sinu_save_folder)

for st in range(len(structure)):
    
    if structure[st] == 'arteriole':
        xl_path = '{}/transformed_arterioles'.format(transformed_structure_folder)
    else:
        xl_path = '{}/transformed_sinosoids'.format(transformed_structure_folder)
        
    xl_data_list = []
    
    for file in os.listdir(xl_path):
        xl_data_list.append(file)
    
    for i in range(len(xl_data_list)):
        print('{} started'.format(xl_data_list[i].split('.')[0]))
            
        if structure[st] == 'arteriole':
            sample_id = xl_data_list[i].split('_A')[0]
        else:
            sample_id = xl_data_list[i].split('_S')[0]
        if structure[st] == 'arteriole':
            if os.path.exists('{}/{}.csv'.format(art_save_folder,sample_id)):
                print('{}_A already exists'.format(sample_id))
                continue
        else:
            if os.path.exists('{}/{}.csv'.format(sinu_save_folder,sample_id)):
                print('{}_S already exists'.format(sample_id))
                continue
        
        obj = sd.read_zarr('{}/{}_with_it_cellularity.zarr'.format(data_folder,sample_id))
        width = int(obj['morphology_focus']['scale0'].dims['x'])
        height = int(obj['morphology_focus']['scale0'].dims['y'])
        del obj
        
        img = Image.new('L', (width, height), 0)
        draw = ImageDraw.Draw(img)
        
        xl = pd.read_csv('{}/{}'.format(xl_path,xl_data_list[i]))
        
        for iii in range(len(xl['geometry'])):
            slctd_poly_temp = wkt.loads(xl['geometry'][iii])
            adjusted_polygon = [(x, y) for x, y in slctd_poly_temp.exterior.coords]
            draw.polygon(adjusted_polygon, outline=1, fill=1)
        
        binary_mask = np.array(img)
        binary_mask_downsampled = cv2.resize((1-binary_mask),(0,0),fx=1/20,fy=1/20,interpolation=cv2.INTER_NEAREST)
        
        distance_map = dm(binary_mask_downsampled)
        
        x_coords, y_coords = np.meshgrid(range(1,distance_map.shape[1]+1), range(1,distance_map.shape[0]+1), indexing='xy')
        
        x_flat = x_coords.flatten()
        y_flat = y_coords.flatten()
        
        df = pd.DataFrame({'x': x_flat, 'y': y_flat, 'value': distance_map.flatten()})
        
        if structure[st] == 'arteriole':
            df.to_csv('{}/{}.csv'.format(art_save_folder,sample_id), index=False)
        else:
            df.to_csv('{}/{}.csv'.format(sinu_save_folder,sample_id), index=False)
        
        print('{} out of {} has been done'.format(i+1,len(xl_data_list)))
