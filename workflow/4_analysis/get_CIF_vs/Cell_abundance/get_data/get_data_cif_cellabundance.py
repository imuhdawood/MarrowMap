import spatialdata as sd
from openslide import OpenSlide
import os
import pandas as pd
import numpy as np
# import geopandas as gpd
from PIL import Image, ImageDraw
import cv2
import tifffile

save_folder = 'path/to/output'
whole_slide_path = 'path/to/wsi'
coord_path = 'path/to/metadata/struct_celltype'
cif_path = 'path/to/cif'
transform_path = 'path/to/landmarks_matrix'
data_folder = 'path/to/main_data'

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

cell_ft_xl = pd.read_csv('../Stats.csv')
fts_grp = cell_ft_xl.columns[4:]
abund_idxs = []
for i in range(len(fts_grp)):
    if not any(char.isdigit() for char in fts_grp[i]):
        abund_idxs.append(fts_grp[i])

abund_idxs.append('CIF')

header = abund_idxs

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    
    mpn_type = cell_ft_xl['diagnosis2'][np.where(np.array(cell_ft_xl['sample_key'])==sample_id)[0][0]]
    
    sub_savepath = '{}/{}'.format(save_folder,mpn_type)
    
    if not os.path.exists(sub_savepath):
        os.makedirs(sub_savepath)
    
    if os.path.exists('{}/{}.csv'.format(sub_savepath,sample_id)):
        print('{} already done'.format(sample_id))
        continue
    
    print('Started {}'.format(sample_id))
    obj = sd.read_zarr(os.path.join(data_folder,'july_sdata_with_adipocytes', zarr_[i]))
    
    its_poly_grp = obj['intertrabecular_regions']['geometry']
    
    trnsmatrix = np.array(pd.read_csv('{}/{}_matrix.csv'.format(transform_path,sample_id),header=None))
    inv_trns = np.linalg.inv(trnsmatrix)
    
    try:
        wsi_img = OpenSlide('{}/{}.ndpi'.format(whole_slide_path,sample_id.split('_')[0]))
    except:
        wsi_img = OpenSlide('{}/{}.ome.tif'.format(whole_slide_path,sample_id.split('_')[0]))
    
    sizes_xy = wsi_img.dimensions
    
    width = sizes_xy[0]
    height = sizes_xy[1]
    
    img = Image.new('L', (width, height), 0)  # 'L' mode for grayscale image

    draw = ImageDraw.Draw(img)
    
    for name,its_poly in its_poly_grp.items():
        lbl_id = int(name.split('_')[-1].split('R')[1])+1
        
        polygons = [(x, y) for x, y in its_poly.exterior.coords]
        adj_polygons = []
        for iii in range(len(polygons)):
            poly_mat = [polygons[iii][0],polygons[iii][1],1]
            trns_x,trns_y,_ = np.dot(inv_trns, poly_mat)
            adj_polygons.append((trns_x,trns_y))
        draw.polygon(adj_polygons, outline=1, fill=lbl_id)
    
    img = np.array(img)
    
    cif_map = tifffile.imread('{}/{}.ome.tif'.format(cif_path,sample_id.split('_')[0]))
    
    img = cv2.resize(img,(cif_map.shape[1],cif_map.shape[0]),interpolation=cv2.INTER_NEAREST)
    
    df_temp = np.zeros((len(its_poly_grp),len(abund_idxs)))
    df_temp = df_temp.astype(float)
    
    for ii in range(len(its_poly_grp)):
        
        its_idx = np.where((np.array(cell_ft_xl['it_regions'])==ii)&(np.array(cell_ft_xl['sample_key'])==sample_id))[0][0]
        cell_its = cell_ft_xl.iloc[its_idx,4:]
        for iii in range(len(abund_idxs)-1):
            df_temp[ii,iii] = cell_its[abund_idxs[iii]]
        lbl_idx = np.where(np.array(img) == (ii+1))
        its_avg_cif = np.mean(cif_map[lbl_idx])
        df_temp[ii,-1] = its_avg_cif
        # pdb.set_trace()
    
    df = pd.DataFrame(df_temp)
    df.columns=header
    
    df.to_csv('{}/{}.csv'.format(sub_savepath,sample_id), index=False)
    
    print('{} has been done'.format(sample_id))
