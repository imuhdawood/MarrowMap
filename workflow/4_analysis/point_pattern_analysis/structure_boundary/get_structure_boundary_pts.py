import spatialdata as sd
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely

save_folder = 'path/to/structure_boundary_points'
tissue_mask_path = 'path/to/tissue_mask'
art_folder = 'path/to/arteriole
sinu_folder = 'path/to/sinusoids'
data_folder = 'path/to/main_data'

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0)

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    if os.path.exists('{}/{}.xlsx'.format(save_folder,sample_id)):
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
    
    bone_strct = obj['transformed_bone']
    
    fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(sample_id) in x]
    fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]
    
    if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
    
        art_gdf_temp = pd.read_csv('{}/{}_A.csv'.format(art_folder,sample_id))
    
        art_gdf = gpd.GeoDataFrame(art_gdf_temp, geometry=shapely.wkt.loads(art_gdf_temp['geometry']))
    
    if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
        sinu_gdf_temp = pd.read_csv('{}/{}_S.csv'.format(sinu_folder,sample_id))
    
        sinu_gdf = gpd.GeoDataFrame(sinu_gdf_temp, geometry=shapely.wkt.loads(sinu_gdf_temp['geometry']))
    
    bone_coord_x = []
    bone_coord_y = []
    
    for bb in range(len(bone_strct)):
        try:
            len_bone = len(bone_strct['geometry'][bb].exterior.coords)
        except:
            geoms = [g for g in bone_strct['geometry'][bb].geoms] #List all multipolygon polygons
            maxx = max(geoms, key=lambda x: x.bounds[2])
            len_bone = len(maxx.exterior.coords)
        pt_step = int(len_bone/10)
        if pt_step < 1:
            pt_step = 1
        try:
            for bbb in range(0,len_bone,pt_step):
                bone_coord_x.append(bone_strct['geometry'][bb].exterior.coords[bbb][0])
                bone_coord_y.append(bone_strct['geometry'][bb].exterior.coords[bbb][1])
        except:
            for bbb in range(0,len_bone,pt_step):
                bone_coord_x.append(maxx.exterior.coords[bbb][0])
                bone_coord_y.append(maxx.exterior.coords[bbb][1])
                
    df_bone = pd.DataFrame(bone_coord_x)
    df_bone.columns =['x']
    df_bone['y'] = bone_coord_y
    
    fat_coord_x = []
    fat_coord_y = []
    
    for bb in range(len(fat_gdf)):
        try:
            len_fat = len(fat_gdf['geometry'][bb].exterior.coords)
        except:
            geoms = [g for g in fat_gdf['geometry'][bb].geoms] #List all multipolygon polygons
            maxx = max(geoms, key=lambda x: x.bounds[2])
            len_fat = len(maxx.exterior.coords)
        pt_step = int(len_fat/10)
        if pt_step < 1:
            pt_step = 1
        try:
            for bbb in range(0,len_fat,pt_step):
                fat_coord_x.append(fat_gdf['geometry'][bb].exterior.coords[bbb][0])
                fat_coord_y.append(fat_gdf['geometry'][bb].exterior.coords[bbb][1])
        except:
            for bbb in range(0,len_fat,pt_step):
                fat_coord_x.append(maxx.exterior.coords[bbb][0])
                fat_coord_y.append(maxx.exterior.coords[bbb][1])
    
    df_fat = pd.DataFrame(fat_coord_x)
    df_fat.columns =['x']
    df_fat['y'] = fat_coord_y
    
    if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
        art_coord_x = []
        art_coord_y = []
        
        for bb in range(len(art_gdf)):
            try:
                len_art = len(art_gdf['geometry'][bb].exterior.coords)
            except:
                geoms = [g for g in art_gdf['geometry'][bb].geoms] #List all multipolygon polygons
                maxx = max(geoms, key=lambda x: x.bounds[2])
                len_art = len(maxx.exterior.coords)
                
            pt_step = int(len_art/10)
            if pt_step < 1:
                pt_step = 1
            try:
                for bbb in range(0,len_art,pt_step):
                    art_coord_x.append(art_gdf['geometry'][bb].exterior.coords[bbb][0])
                    art_coord_y.append(art_gdf['geometry'][bb].exterior.coords[bbb][1])
            except:
                for bbb in range(0,len_art,pt_step):
                    art_coord_x.append(maxx.exterior.coords[bbb][0])
                    art_coord_y.append(maxx.exterior.coords[bbb][1])
        
        df_art = pd.DataFrame(art_coord_x)
        df_art.columns =['x']
        df_art['y'] = art_coord_y
    else:
        art_coord_x = [np.nan]
        art_coord_y = [np.nan]
        df_art = pd.DataFrame(art_coord_x)
        df_art.columns =['x']
        df_art['y'] = art_coord_y
    
    if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
        sinu_coord_x = []
        sinu_coord_y = []
        
        for bb in range(len(sinu_gdf)):
            try:
                len_sinu = len(sinu_gdf['geometry'][bb].exterior.coords)
            except:
                geoms = [g for g in sinu_gdf['geometry'][bb].geoms] #List all multipolygon polygons
                maxx = max(geoms, key=lambda x: x.bounds[2])
                len_sinu = len(maxx.exterior.coords)
            pt_step = int(len_sinu/10)
            if pt_step < 1:
                pt_step = 1
            try:
                for bbb in range(0,len_sinu,pt_step):
                    sinu_coord_x.append(sinu_gdf['geometry'][bb].exterior.coords[bbb][0])
                    sinu_coord_y.append(sinu_gdf['geometry'][bb].exterior.coords[bbb][1])
            except:
                geoms = [g for g in sinu_gdf['geometry'][bb].geoms] #List all multipolygon polygons
                maxx = max(geoms, key=lambda x: x.bounds[2])
                for bbb in range(0,len_sinu,pt_step):
                    sinu_coord_x.append(maxx.exterior.coords[bbb][0])
                    sinu_coord_y.append(maxx.exterior.coords[bbb][1])
        
        df_sinu = pd.DataFrame(sinu_coord_x)
        df_sinu.columns =['x']
        df_sinu['y'] = sinu_coord_y
    else:
        sinu_coord_x = [np.nan]
        sinu_coord_y = [np.nan]
        df_sinu = pd.DataFrame(sinu_coord_x)
        df_sinu.columns =['x']
        df_sinu['y'] = sinu_coord_y
    
    writer_strct = pd.ExcelWriter('{}/{}.xlsx'.format(save_folder,sample_id), engine='xlsxwriter')
    df_bone.to_excel(writer_strct, sheet_name="Bone", index=False)
    df_fat.to_excel(writer_strct, sheet_name='Fat', index=False)
    df_art.to_excel(writer_strct, sheet_name='Arteriole', index=False)
    df_sinu.to_excel(writer_strct, sheet_name='Sinusoid', index=False)
    
    writer_strct.close()
    
    print('{} out of {} : {} has been done'.format(i+1,len(zarr_),sample_id))
