import spatialdata as sd
import os
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import box
import shapely

save_folder = 'path/to/output'
art_folder = 'path/to/arteriole' # Transformed arteriole
sinu_folder = 'path/to/arteriole' # Transformed sinusoid
data_folder = 'path/to/main_data'
dim_folder = 'path/to/dimension'

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

zarr_ = []
for file in os.listdir('{}/july_sdata_with_adipocytes'.format(data_folder)):
    if file.endswith('.zarr'):
        zarr_.append(file)

df = pd.read_csv('{}/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv'.format(data_folder), index_col=0)

def filter_within_bounds(gdf, min_x, min_y, max_x, max_y):
    bbox = box(min_x, min_y, max_x, max_y)
    return gdf[gdf.intersects(bbox)]

def find_next_nearest(row, gdf):
    poly_id = row["index"]
    poly_geom = row.geometry

    # Exclude self and compute distances
    gdf_other = gdf[gdf["index"] != poly_id].copy()
    gdf_other["distance"] = gdf_other.geometry.distance(poly_geom)
    
    # Get the closest valid neighbor
    nearest_row = gdf_other.nsmallest(1, "distance")
    
    if nearest_row.empty:
        return None  # Handle case where no neighbor is found
    
    return {
        "geometry": poly_geom,  # Original polygon geometry
        "nearest_poly_id": nearest_row["index"].values[0],  # Nearest neighbor's ID
        "nearest_geometry": nearest_row.geometry.values[0],  # Nearest neighbor's geometry
        "distance": nearest_row["distance"].values[0],  # Distance to nearest neighbor
    }

for i in range(len(zarr_)):
    sample_id = zarr_[i].split('_with')[0]
    
    if not os.path.exists('{}/{}.xlsx'.format(dim_folder,sample_id)):
        print('{} out of {} : {} has been done'.format(i+1,len(zarr_),sample_id))
        continue
    
    if os.path.exists('{}/{}.csv'.format(save_folder,sample_id)):
        print('{} already done'.format(sample_id))
        continue
    
    dim_xl = pd.read_excel('{}/{}.xlsx'.format(dim_folder,sample_id))
    
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
    
    art_gdf_temp = pd.read_csv('{}/{}_A.csv'.format(art_folder,sample_id))

    art_gdf = gpd.GeoDataFrame(art_gdf_temp, geometry=shapely.wkt.loads(art_gdf_temp['geometry']))
    
    sinu_gdf_temp = pd.read_csv('{}/{}_S.csv'.format(sinu_folder,sample_id))

    sinu_gdf = gpd.GeoDataFrame(sinu_gdf_temp, geometry=shapely.wkt.loads(sinu_gdf_temp['geometry']))
    
    bone_gdf_filtered = filter_within_bounds(bone_strct,dim_xl['min_x'],dim_xl['min_y'],dim_xl['max_x'],dim_xl['max_y'])
    
    fat_gdf_filtered = filter_within_bounds(fat_gdf,dim_xl['min_x'],dim_xl['min_y'],dim_xl['max_x'],dim_xl['max_y'])
    
    art_gdf_filtered = filter_within_bounds(art_gdf,dim_xl['min_x'],dim_xl['min_y'],dim_xl['max_x'],dim_xl['max_y'])
    
    sinu_gdf_filtered = filter_within_bounds(sinu_gdf,dim_xl['min_x'],dim_xl['min_y'],dim_xl['max_x'],dim_xl['max_y'])
    
    bone_gdf_filtered["index"] = bone_gdf_filtered.index
    fat_gdf_filtered["index"] = fat_gdf_filtered.index 
    art_gdf_filtered["index"] = art_gdf_filtered.index 
    sinu_gdf_filtered["index"] = sinu_gdf_filtered.index 
    
    distance_bone_to_bone = []
    
    if len(bone_gdf_filtered) > 1:
        bone_to_bone_temp = [find_next_nearest(row, bone_gdf_filtered) for _, row in bone_gdf_filtered.iterrows()]
        bone_to_bone_temp = [res for res in bone_to_bone_temp if res is not None]  # Remove None values
        
        # Convert to GeoDataFrame
        bone_to_bone = gpd.GeoDataFrame(bone_to_bone_temp, geometry="geometry", index=bone_gdf_filtered["index"])
        
        bone_to_bone_v = bone_to_bone.rename(columns={"nearest_poly_id": "bone_to_nearest_bone", "distance": "distance_to_bone"})
        
        bone_to_bone_v = bone_to_bone_v[["bone_to_nearest_bone", "distance_to_bone"]]
        
        for name, poly in bone_to_bone_v['distance_to_bone'].items():
            distance_bone_to_bone.append(poly)
    else:
        distance_bone_to_bone.append(0)
    
    bone_to_fat = gpd.sjoin_nearest(bone_gdf_filtered, fat_gdf_filtered, distance_col="distances", how="left")
    bone_to_fat_dropped = bone_to_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    bone_to_fat_dropped = bone_to_fat_dropped.rename(columns={"index_right": "bone_to_nearest_fat", "distances": "distance_to_fat"})
    bone_to_fat_dropped = bone_to_fat_dropped[["bone_to_nearest_fat", "distance_to_fat"]]
    
    distance_bone_to_fat = []
    for name, poly in bone_to_fat_dropped['distance_to_fat'].items():
        distance_bone_to_fat.append(poly)
    
    distance_bone_to_art = []
    
    bone_to_art = gpd.sjoin_nearest(bone_gdf_filtered, art_gdf_filtered, distance_col="distances", how="left")
    bone_to_art_dropped = bone_to_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    bone_to_art_dropped = bone_to_art_dropped.rename(columns={"index_right": "bone_to_nearest_art", "distances": "distance_to_art"})
    bone_to_art_dropped = bone_to_art_dropped[["bone_to_nearest_art", "distance_to_art"]]

    # distance_bone_to_art = []
    for name, poly in bone_to_art_dropped['distance_to_art'].items():
        distance_bone_to_art.append(poly)
    
    distance_bone_to_sinu = []
    
    bone_to_sinu = gpd.sjoin_nearest(bone_gdf_filtered, sinu_gdf_filtered, distance_col="distances", how="left")
    bone_to_sinu_dropped = bone_to_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    bone_to_sinu_dropped = bone_to_sinu_dropped.rename(columns={"index_right": "bone_to_nearest_sinu", "distances": "distance_to_sinu"})
    bone_to_sinu_dropped = bone_to_sinu_dropped[["bone_to_nearest_sinu", "distance_to_sinu"]]

    # distance_bone_to_sinu = []
    for name, poly in bone_to_sinu_dropped['distance_to_sinu'].items():
        distance_bone_to_sinu.append(poly)
    
    fat_to_bone = gpd.sjoin_nearest(fat_gdf_filtered, bone_gdf_filtered, distance_col="distances", how="left")
    fat_to_bone_dropped = fat_to_bone.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    fat_to_bone_dropped = fat_to_bone_dropped.rename(columns={"index_right": "fat_to_nearest_bone", "distances": "distance_to_bone"})
    fat_to_bone_dropped = fat_to_bone_dropped[["fat_to_nearest_bone", "distance_to_bone"]]
    
    distance_fat_to_bone = []
    for name, poly in fat_to_bone_dropped['distance_to_bone'].items():
        distance_fat_to_bone.append(poly)
    
    distance_fat_to_fat = []
    if len(fat_gdf_filtered) > 1:
        fat_to_fat_temp = [find_next_nearest(row, fat_gdf_filtered) for _, row in fat_gdf_filtered.iterrows()]
        fat_to_fat_temp = [res for res in fat_to_fat_temp if res is not None]  # Remove None values
        
        # Convert to GeoDataFrame
        fat_to_fat = gpd.GeoDataFrame(fat_to_fat_temp, geometry="geometry", index=fat_gdf_filtered["index"])
        
        fat_to_fat_v = fat_to_fat.rename(columns={"nearest_poly_id": "fat_to_nearest_fat", "distance": "distance_to_fat"})
        
        fat_to_fat_v = fat_to_fat_v[["fat_to_nearest_fat", "distance_to_fat"]]
        
        for name, poly in fat_to_fat_v['distance_to_fat'].items():
            distance_fat_to_fat.append(poly)
    else:
        distance_fat_to_fat.append(0)
        
    distance_fat_to_art = []
    
    fat_to_art = gpd.sjoin_nearest(fat_gdf_filtered, art_gdf_filtered, distance_col="distances", how="left")
    fat_to_art_dropped = fat_to_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    fat_to_art_dropped = fat_to_art_dropped.rename(columns={"index_right": "fat_to_nearest_art", "distances": "distance_to_art"})
    fat_to_art_dropped = fat_to_art_dropped[["fat_to_nearest_art", "distance_to_art"]]

    # distance_fat_to_art = []
    for name, poly in fat_to_art_dropped['distance_to_art'].items():
        distance_fat_to_art.append(poly)
    
    distance_fat_to_sinu = []
    
    fat_to_sinu = gpd.sjoin_nearest(fat_gdf_filtered, sinu_gdf_filtered, distance_col="distances", how="left")
    fat_to_sinu_dropped = fat_to_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    fat_to_sinu_dropped = fat_to_sinu_dropped.rename(columns={"index_right": "fat_to_nearest_sinu", "distances": "distance_to_sinu"})
    fat_to_sinu_dropped = fat_to_sinu_dropped[["fat_to_nearest_sinu", "distance_to_sinu"]]

    # distance_fat_to_sinu = []
    for name, poly in fat_to_sinu_dropped['distance_to_sinu'].items():
        distance_fat_to_sinu.append(poly)
        
    distance_art_to_bone = []
    distance_art_to_fat = []
    distance_art_to_sinu = []
    distance_art_to_art = []
    
    art_to_bone = gpd.sjoin_nearest(art_gdf_filtered, bone_gdf_filtered, distance_col="distances", how="left")
    art_to_bone_dropped = art_to_bone.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    art_to_bone_dropped = art_to_bone_dropped.rename(columns={"index_right": "art_to_nearest_bone", "distances": "distance_to_bone"})
    art_to_bone_dropped = art_to_bone_dropped[["art_to_nearest_bone", "distance_to_bone"]]
    
    # distance_art_to_bone = []
    for name, poly in art_to_bone_dropped['distance_to_bone'].items():
        distance_art_to_bone.append(poly)
    
    art_to_fat = gpd.sjoin_nearest(art_gdf_filtered, fat_gdf_filtered, distance_col="distances", how="left")
    art_to_fat_dropped = art_to_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    art_to_fat_dropped = art_to_fat_dropped.rename(columns={"index_right": "art_to_nearest_fat", "distances": "distance_to_fat"})
    art_to_fat_dropped = art_to_fat_dropped[["art_to_nearest_fat", "distance_to_fat"]]
    
    # distance_art_to_fat = []
    for name, poly in art_to_fat_dropped['distance_to_fat'].items():
        distance_art_to_fat.append(poly)
    
    distance_art_to_art = []
    if len(art_gdf_filtered) > 1:
        art_to_art_temp = [find_next_nearest(row, art_gdf_filtered) for _, row in art_gdf_filtered.iterrows()]
        art_to_art_temp = [res for res in art_to_art_temp if res is not None]  # Remove None values
        
        # Convert to GeoDataFrame
        art_to_art = gpd.GeoDataFrame(art_to_art_temp, geometry="geometry", index=art_gdf_filtered["index"])
        
        art_to_art_v = art_to_art.rename(columns={"nearest_poly_id": "art_to_nearest_art", "distance": "distance_to_art"})
        
        art_to_art_v = art_to_art_v[["art_to_nearest_art", "distance_to_art"]]
        
        for name, poly in art_to_art_v['distance_to_art'].items():
            distance_art_to_art.append(poly)
    else:
        distance_art_to_art.append(0)
    
    art_to_sinu = gpd.sjoin_nearest(art_gdf_filtered, sinu_gdf_filtered, distance_col="distances", how="left")
    art_to_sinu_dropped = art_to_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    art_to_sinu_dropped = art_to_sinu_dropped.rename(columns={"index_right": "art_to_nearest_sinu", "distances": "distance_to_sinu"})
    art_to_sinu_dropped = art_to_sinu_dropped[["art_to_nearest_sinu", "distance_to_sinu"]]

    # distance_art_to_sinu = []
    for name, poly in art_to_sinu_dropped['distance_to_sinu'].items():
        distance_art_to_sinu.append(poly)
    
    distance_sinu_to_bone = []
    distance_sinu_to_fat = []
    distance_sinu_to_art = []
    distance_sinu_to_sinu = []
    
    sinu_to_bone = gpd.sjoin_nearest(sinu_gdf_filtered, bone_gdf_filtered, distance_col="distances", how="left")
    sinu_to_bone_dropped = sinu_to_bone.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    sinu_to_bone_dropped = sinu_to_bone_dropped.rename(columns={"index_right": "sinu_to_nearest_bone", "distances": "distance_to_bone"})
    sinu_to_bone_dropped = sinu_to_bone_dropped[["sinu_to_nearest_bone", "distance_to_bone"]]
    
    for name, poly in sinu_to_bone_dropped['distance_to_bone'].items():
        distance_sinu_to_bone.append(poly)
    
    sinu_to_fat = gpd.sjoin_nearest(sinu_gdf_filtered, fat_gdf_filtered, distance_col="distances", how="left")
    sinu_to_fat_dropped = sinu_to_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    sinu_to_fat_dropped = sinu_to_fat_dropped.rename(columns={"index_right": "sinu_to_nearest_fat", "distances": "distance_to_fat"})
    sinu_to_fat_dropped = sinu_to_fat_dropped[["sinu_to_nearest_fat", "distance_to_fat"]]
    
    for name, poly in sinu_to_fat_dropped['distance_to_fat'].items():
        distance_sinu_to_fat.append(poly)
    
    sinu_to_art = gpd.sjoin_nearest(sinu_gdf_filtered, art_gdf_filtered, distance_col="distances", how="left")
    sinu_to_art_dropped = sinu_to_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
    sinu_to_art_dropped = sinu_to_art_dropped.rename(columns={"index_right": "sinu_to_nearest_art", "distances": "distance_to_art"})
    sinu_to_art_dropped = sinu_to_art_dropped[["sinu_to_nearest_art", "distance_to_art"]]

    for name, poly in sinu_to_art_dropped['distance_to_art'].items():
        distance_sinu_to_art.append(poly)
    
    distance_sinu_to_sinu = []
    if len(sinu_gdf_filtered) > 1:
        sinu_to_sinu_temp = [find_next_nearest(row, sinu_gdf_filtered) for _, row in sinu_gdf_filtered.iterrows()]
        sinu_to_sinu_temp = [res for res in sinu_to_sinu_temp if res is not None]  # Remove None values
        
        # Convert to GeoDataFrame
        sinu_to_sinu = gpd.GeoDataFrame(sinu_to_sinu_temp, geometry="geometry", index=sinu_gdf_filtered["index"])
        
        sinu_to_sinu_v = sinu_to_sinu.rename(columns={"nearest_poly_id": "sinu_to_nearest_sinu", "distance": "distance_to_sinu"})
        
        sinu_to_sinu_v = sinu_to_sinu_v[["sinu_to_nearest_sinu", "distance_to_sinu"]]
        
        distance_sinu_to_sinu = []
        for name, poly in sinu_to_sinu_v['distance_to_sinu'].items():
            distance_sinu_to_sinu.append(poly)
    else:
        distance_sinu_to_sinu.append(0)
    
    max_lenth = max(len(distance_bone_to_fat),len(distance_bone_to_art),len(distance_bone_to_sinu),len(distance_fat_to_bone),len(distance_fat_to_art),len(distance_fat_to_sinu),len(distance_art_to_bone),len(distance_art_to_fat),len(distance_art_to_sinu),len(distance_sinu_to_bone),len(distance_sinu_to_fat),len(distance_sinu_to_art))
    
    header = ['Bone to Bone','Bone to Fat', 'Bone to Arteriole', 'Bone to Sinusoid', 
              'Fat to Bone', 'Fat to Fat', 'Fat to Arteriole', 'Fat to Sinusoid', 
              'Arteriole to Bone', 'Arteriole to Fat', 'Arteriole to Arteriole', 'Arteriole to Sinusoid', 
              'Sinusoid to Bone', 'Sinusoid to Fat', 'Sinusoid to Arteriole', 'Sinusoid to Sinusoid']
    
    dff = pd.DataFrame(np.nan, index=range(max_lenth), columns=header)
    dff.loc[:len(distance_bone_to_bone)-1, 'Bone to Bone'] = distance_bone_to_bone
    dff.loc[:len(distance_bone_to_fat)-1, 'Bone to Fat'] = distance_bone_to_fat
    dff.loc[:len(distance_bone_to_art)-1, 'Bone to Arteriole'] = distance_bone_to_art
    dff.loc[:len(distance_bone_to_sinu)-1, 'Bone to Sinusoid'] = distance_bone_to_sinu
    dff.loc[:len(distance_fat_to_bone)-1, 'Fat to Bone'] = distance_fat_to_bone
    dff.loc[:len(distance_fat_to_fat)-1, 'Fat to Fat'] = distance_fat_to_fat
    dff.loc[:len(distance_fat_to_art)-1, 'Fat to Arteriole'] = distance_fat_to_art
    dff.loc[:len(distance_fat_to_sinu)-1, 'Fat to Sinusoid'] = distance_fat_to_sinu
    dff.loc[:len(distance_art_to_bone)-1, 'Arteriole to Bone'] = distance_art_to_bone
    dff.loc[:len(distance_art_to_fat)-1, 'Arteriole to Fat'] = distance_art_to_fat
    dff.loc[:len(distance_art_to_art)-1, 'Arteriole to Arteriole'] = distance_art_to_art
    dff.loc[:len(distance_art_to_sinu)-1, 'Arteriole to Sinusoid'] = distance_art_to_sinu
    dff.loc[:len(distance_sinu_to_bone)-1, 'Sinusoid to Bone'] = distance_sinu_to_bone
    dff.loc[:len(distance_sinu_to_fat)-1, 'Sinusoid to Fat'] = distance_sinu_to_fat
    dff.loc[:len(distance_sinu_to_art)-1, 'Sinusoid to Arteriole'] = distance_sinu_to_art
    dff.loc[:len(distance_sinu_to_sinu)-1, 'Sinusoid to Sinusoid'] = distance_sinu_to_sinu
    
    dff.to_csv('{}/{}.csv'.format(save_folder,sample_id), index=False)
    print('{} out of {} : {} has been done'.format(i+1,len(zarr_),sample_id))
