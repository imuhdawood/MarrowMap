import os
import spatialdata as sd
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely
from shapely.geometry import box

data_folder = 'path/to/data'
main_savepth = 'metadata'

compare_type = ['struct_struct', 'struct_celltype', 'struct_cellneighbor',
                'celltype_celltype','cellneighbor_cellneighbor']

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
    
    for comptype in compare_type:
        if 'struct' in comptype:
            
            art_folder = 'path/to/arteriole' # Transformed arteriole
            sinu_folder = 'path/to/sinusoid' # Transformed sinusoid
            
            bone_strct = obj['transformed_bone']
            
            fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(sample_id) in x]
            fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]
            
            art_gdf_temp = pd.read_csv('{}/{}_A.csv'.format(art_folder,sample_id))
            art_gdf = gpd.GeoDataFrame(art_gdf_temp, geometry=shapely.wkt.loads(art_gdf_temp['geometry']))
            
            sinu_gdf_temp = pd.read_csv('{}/{}_S.csv'.format(sinu_folder,sample_id))
            sinu_gdf = gpd.GeoDataFrame(sinu_gdf_temp, geometry=shapely.wkt.loads(sinu_gdf_temp['geometry']))
            
            if comptype.split['_'][1] == 'struct':
                
                dim_folder = 'path/to/dimension'
                save_folder = '{}/{}'.format(main_savepth,comptype)
                
                if os.path.exists('{}/{}.xlsx'.format(dim_folder,sample_id)):
                    dim_xl = pd.read_excel('{}/{}.xlsx'.format(dim_folder,sample_id))
                else:
                    dim_xl = pd.DataFrame([{
                                            "min_x": 0,
                                            "min_y": 0,
                                            "max_x": int(obj['morphology_focus']['scale0'].dims['x']),
                                            "max_y": int(obj['morphology_focus']['scale0'].dims['y'])
                                        }])
                
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
                
                
                distance_bone_to_art = []
                distance_fat_to_art = []
                distance_art_to_bone = []
                distance_art_to_fat = []
                distance_art_to_sinu = []
                distance_art_to_art = []
                
                if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
                
                    bone_to_art = gpd.sjoin_nearest(bone_gdf_filtered, art_gdf_filtered, distance_col="distances", how="left")
                    bone_to_art_dropped = bone_to_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    bone_to_art_dropped = bone_to_art_dropped.rename(columns={"index_right": "bone_to_nearest_art", "distances": "distance_to_art"})
                    bone_to_art_dropped = bone_to_art_dropped[["bone_to_nearest_art", "distance_to_art"]]
    
                    for name, poly in bone_to_art_dropped['distance_to_art'].items():
                        distance_bone_to_art.append(poly)
                    
                    fat_to_art = gpd.sjoin_nearest(fat_gdf_filtered, art_gdf_filtered, distance_col="distances", how="left")
                    fat_to_art_dropped = fat_to_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    fat_to_art_dropped = fat_to_art_dropped.rename(columns={"index_right": "fat_to_nearest_art", "distances": "distance_to_art"})
                    fat_to_art_dropped = fat_to_art_dropped[["fat_to_nearest_art", "distance_to_art"]]

                    for name, poly in fat_to_art_dropped['distance_to_art'].items():
                        distance_fat_to_art.append(poly)
                    
                    art_to_bone = gpd.sjoin_nearest(art_gdf_filtered, bone_gdf_filtered, distance_col="distances", how="left")
                    art_to_bone_dropped = art_to_bone.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    art_to_bone_dropped = art_to_bone_dropped.rename(columns={"index_right": "art_to_nearest_bone", "distances": "distance_to_bone"})
                    art_to_bone_dropped = art_to_bone_dropped[["art_to_nearest_bone", "distance_to_bone"]]
                    
                    for name, poly in art_to_bone_dropped['distance_to_bone'].items():
                        distance_art_to_bone.append(poly)
                    
                    art_to_fat = gpd.sjoin_nearest(art_gdf_filtered, fat_gdf_filtered, distance_col="distances", how="left")
                    art_to_fat_dropped = art_to_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    art_to_fat_dropped = art_to_fat_dropped.rename(columns={"index_right": "art_to_nearest_fat", "distances": "distance_to_fat"})
                    art_to_fat_dropped = art_to_fat_dropped[["art_to_nearest_fat", "distance_to_fat"]]
                    
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
                    
                    if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
                    
                        art_to_sinu = gpd.sjoin_nearest(art_gdf_filtered, sinu_gdf_filtered, distance_col="distances", how="left")
                        art_to_sinu_dropped = art_to_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                        art_to_sinu_dropped = art_to_sinu_dropped.rename(columns={"index_right": "art_to_nearest_sinu", "distances": "distance_to_sinu"})
                        art_to_sinu_dropped = art_to_sinu_dropped[["art_to_nearest_sinu", "distance_to_sinu"]]
                        
                        for name, poly in art_to_sinu_dropped['distance_to_sinu'].items():
                            distance_art_to_sinu.append(poly)
                    
                
                distance_bone_to_sinu = []
                distance_fat_to_sinu = []
                distance_sinu_to_bone = []
                distance_sinu_to_fat = []
                distance_sinu_to_art = []
                distance_sinu_to_sinu = []
                
                if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
                    distance_bone_to_sinu = []
                    
                    bone_to_sinu = gpd.sjoin_nearest(bone_gdf_filtered, sinu_gdf_filtered, distance_col="distances", how="left")
                    bone_to_sinu_dropped = bone_to_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    bone_to_sinu_dropped = bone_to_sinu_dropped.rename(columns={"index_right": "bone_to_nearest_sinu", "distances": "distance_to_sinu"})
                    bone_to_sinu_dropped = bone_to_sinu_dropped[["bone_to_nearest_sinu", "distance_to_sinu"]]

                    for name, poly in bone_to_sinu_dropped['distance_to_sinu'].items():
                        distance_bone_to_sinu.append(poly)
                    
                    fat_to_sinu = gpd.sjoin_nearest(fat_gdf_filtered, sinu_gdf_filtered, distance_col="distances", how="left")
                    fat_to_sinu_dropped = fat_to_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    fat_to_sinu_dropped = fat_to_sinu_dropped.rename(columns={"index_right": "fat_to_nearest_sinu", "distances": "distance_to_sinu"})
                    fat_to_sinu_dropped = fat_to_sinu_dropped[["fat_to_nearest_sinu", "distance_to_sinu"]]

                    for name, poly in fat_to_sinu_dropped['distance_to_sinu'].items():
                        distance_fat_to_sinu.append(poly)
                    
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
                    
                    if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
                        sinu_to_art = gpd.sjoin_nearest(sinu_gdf_filtered, art_gdf_filtered, distance_col="distances", how="left")
                        sinu_to_art_dropped = sinu_to_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                        sinu_to_art_dropped = sinu_to_art_dropped.rename(columns={"index_right": "sinu_to_nearest_art", "distances": "distance_to_art"})
                        sinu_to_art_dropped = sinu_to_art_dropped[["sinu_to_nearest_art", "distance_to_art"]]
    
                        for name, poly in sinu_to_art_dropped['distance_to_art'].items():
                            distance_sinu_to_art.append(poly)
                    
                
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
                print('{}: {} out of {} : {} has been done'.format(comptype, i+1,len(zarr_),sample_id))
                
                del dff
                del distance_bone_to_bone
                del distance_bone_to_fat
                del distance_bone_to_art
                del distance_bone_to_sinu
                del distance_fat_to_bone
                del distance_fat_to_fat
                del distance_fat_to_art
                del distance_fat_to_sinu
                del distance_art_to_bone
                del distance_art_to_fat
                del distance_art_to_art
                del distance_art_to_sinu
                del distance_sinu_to_bone
                del distance_sinu_to_fat
                del distance_sinu_to_art
                del distance_sinu_to_sinu
                
                
            elif comptype.split['_'][1] == 'celltype':
                
                save_folder = '{}/{}'.format(main_savepth,comptype)
                
                cell_centroid = obj["new_cell_boundaries"]['geometry'].centroid
                cell_ids = adata.obs["cell_id_x"]
                cell_types = adata.obs["obj.anno_3_w_megs_w_stromal_x"]
                cell_bound = obj["new_cell_boundaries"]
                
                cell_centroid_x = []
                cell_centroid_y = []
                
                for name, poly in cell_centroid.items():
                    if name in cell_ids.values:
                        cell_centroid_x.append(poly.x)
                        cell_centroid_y.append(poly.y)
                
                distance_bone = adata.obs['distance_to_bone']
                fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(sample_id) in x]
                fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]
                
                cells_w_fat = gpd.sjoin_nearest(cell_bound, fat_gdf, distance_col="distances", how="left")
                fat_dropped = cells_w_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                fat_dropped = fat_dropped.rename(columns={"index_right": "nearest_fat", "distances": "distance_to_fat"})
                fat_dropped = fat_dropped[["nearest_fat", "distance_to_fat"]]
                
                distance_fat = []
                for name, poly in fat_dropped['distance_to_fat'].items():
                    if name in cell_ids.values:
                        distance_fat.append(poly)
                
                distance_art = []
                
                if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
                
                    art_gdf_temp = pd.read_csv('{}/{}_A.csv'.format(art_folder,sample_id))
                
                    art_gdf = gpd.GeoDataFrame(art_gdf_temp, geometry=shapely.wkt.loads(art_gdf_temp['geometry']))
                    
                    cells_w_art = gpd.sjoin_nearest(cell_bound, art_gdf, distance_col="distances", how="left")
                    art_dropped = cells_w_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    art_dropped = art_dropped.rename(columns={"index_right": "nearest_arteriole", "distances": "distance_to_arteriole"})
                    art_dropped = art_dropped[["nearest_arteriole", "distance_to_arteriole"]]
                    
                    for name, poly in art_dropped['distance_to_arteriole'].items():
                        if name in cell_ids.values:
                            distance_art.append(poly)
                
                distance_sinu = []
                
                if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
                    sinu_gdf_temp = pd.read_csv('{}/{}_S.csv'.format(sinu_folder,sample_id))
                
                    sinu_gdf = gpd.GeoDataFrame(sinu_gdf_temp, geometry=shapely.wkt.loads(sinu_gdf_temp['geometry']))
                    
                    cells_w_sinu = gpd.sjoin_nearest(cell_bound, sinu_gdf, distance_col="distances", how="left")
                    sinu_dropped = cells_w_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    sinu_dropped = sinu_dropped.rename(columns={"index_right": "nearest_sinusoid", "distances": "distance_to_sinusoid"})
                    sinu_dropped = sinu_dropped[["nearest_sinusoid", "distance_to_sinusoid"]]
                    
                    for name, poly in sinu_dropped['distance_to_sinusoid'].items():
                        if name in cell_ids.values:
                            distance_sinu.append(poly)
                
                max_lenth = max(len(cell_types),len(distance_bone),len(distance_fat),len(distance_art),len(distance_sinu))
                
                header = ['Cell Type', 'x', 'y', 'Distance to Bone', 
                          'Distance to Fat', 'Distance to Arteriole', 'Distance to Sinusoid']
                
                dff = pd.DataFrame(np.nan, index=range(max_lenth), columns=header)
                dff.loc[:len(cell_types)-1, 'Cell Type'] = cell_types
                dff.loc[:len(cell_centroid_x)-1, 'x'] = cell_centroid_x
                dff.loc[:len(cell_centroid_y)-1, 'y'] = cell_centroid_y
                dff.loc[:len(distance_bone)-1, 'Distance to Bone'] = distance_bone
                dff.loc[:len(distance_fat)-1, 'Distance to Fat'] = distance_fat
                dff.loc[:len(distance_art)-1, 'Distance to Arteriole'] = distance_art
                dff.loc[:len(distance_sinu)-1, 'Distance to Sinusoid'] = distance_sinu
                dff.to_csv('{}/{}.csv'.format(save_folder,sample_id), index=False)
                
                
                print('{}: {} out of {} : {} has been done'.format(comptype, i+1,len(zarr_),sample_id))
                
                del cell_types
                del cell_centroid_x
                del cell_centroid_y
                del distance_bone
                del distance_fat
                del distance_art
                del distance_sinu
                
            elif comptype.split['_'][1] == 'cellneighbor':
                
                save_folder = '{}/{}'.format(main_savepth,comptype)
                
                cell_centroid = obj["new_cell_boundaries"]['geometry'].centroid
                cell_ids_orig = adata.obs["cell_id_x"]
                cluster_df = pd.read_csv('anndata_cc_3_n10.csv', index_col=0)
                cluster_ids = np.array(cluster_df.index)
                
                cluster_type = []
                cell_ids = []
                
                set_cluster_ids = set(cluster_ids)
                
                in_mask = [name in set_cluster_ids for name in cell_ids_orig]
                
                cell_ids = np.array(cell_ids_orig)[in_mask]
                
                distance_bone = adata.obs['distance_to_bone'].iloc[in_mask]
                
                set_cell_ids = set(cell_ids)
                in_mask = [name in set_cell_ids for name in np.array(cell_centroid.index)]
                cell_centroid_ = cell_centroid[in_mask]
                
                del cell_centroid
                
                cluster_idx = np.where([name in set_cell_ids for name in cluster_ids])[0]
                cluster_type = cluster_df['cluster_cellcharter_n5'].iloc[cluster_idx]
                cell_bound = obj["new_cell_boundaries"]
                
                
                fat_cells = [x for x in obj["new_cell_boundaries_with_full_adipocytes"].index.to_list() if '{}_A'.format(sample_id) in x]
                fat_gdf = obj["new_cell_boundaries_with_full_adipocytes"].loc[fat_cells]
                
                cells_w_fat = gpd.sjoin_nearest(cell_bound, fat_gdf, distance_col="distances", how="left")
                fat_dropped = cells_w_fat.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                fat_dropped = fat_dropped.rename(columns={"index_right": "nearest_fat", "distances": "distance_to_fat"})
                fat_dropped = fat_dropped[["nearest_fat", "distance_to_fat"]]
                
                
                distance_fat = fat_dropped['distance_to_fat'][in_mask]
                
                if os.path.exists('{}/{}_A.csv'.format(art_folder,sample_id)):
                    art_gdf_temp = pd.read_csv('{}/{}_A.csv'.format(art_folder,sample_id))
                
                    art_gdf = gpd.GeoDataFrame(art_gdf_temp, geometry=shapely.wkt.loads(art_gdf_temp['geometry']))
                    
                    cells_w_art = gpd.sjoin_nearest(cell_bound, art_gdf, distance_col="distances", how="left")
                    art_dropped = cells_w_art.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    art_dropped = art_dropped.rename(columns={"index_right": "nearest_arteriole", "distances": "distance_to_arteriole"})
                    art_dropped = art_dropped[["nearest_arteriole", "distance_to_arteriole"]]
                    
                    # in_mask = [name in set_cell_ids for name in np.array(art_dropped['distance_to_arteriole'].index)]
                    
                    distance_art = art_dropped['distance_to_arteriole'][in_mask]
                
                else:
                    distance_art = []
                
                if os.path.exists('{}/{}_S.csv'.format(sinu_folder,sample_id)):
                    sinu_gdf_temp = pd.read_csv('{}/{}_S.csv'.format(sinu_folder,sample_id))
                
                    sinu_gdf = gpd.GeoDataFrame(sinu_gdf_temp, geometry=shapely.wkt.loads(sinu_gdf_temp['geometry']))
                    
                    cells_w_sinu = gpd.sjoin_nearest(cell_bound, sinu_gdf, distance_col="distances", how="left")
                    sinu_dropped = cells_w_sinu.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
                    sinu_dropped = sinu_dropped.rename(columns={"index_right": "nearest_sinusoid", "distances": "distance_to_sinusoid"})
                    sinu_dropped = sinu_dropped[["nearest_sinusoid", "distance_to_sinusoid"]]
                    
                    # in_mask = [name in set_cell_ids for name in np.array(sinu_dropped['distance_to_sinusoid'].index)]
                    
                    distance_sinu = sinu_dropped['distance_to_sinusoid'][in_mask]
                
                else:
                    distance_sinu = []
                
                max_lenth = max(len(cluster_type),len(distance_bone),len(distance_fat),len(distance_art),len(distance_sinu))
                
                header = ['Cluster Type', 'x', 'y', 'Distance to Bone', 'Distance to Fat', 
                          'Distance to Arteriole', 'Distance to Sinusoid']
                
                dff = pd.DataFrame(np.nan, index=range(max_lenth), columns=header)
                dff.loc[:len(cluster_type)-1, 'Cluster Type'] = cluster_type
                dff.loc[:len(cell_centroid_x)-1, 'x'] = np.array(cell_centroid_.x)
                dff.loc[:len(cell_centroid_y)-1, 'y'] = np.array(cell_centroid_.y)
                dff.loc[:len(distance_bone)-1, 'Distance to Bone'] = distance_bone
                dff.loc[:len(distance_fat)-1, 'Distance to Fat'] = distance_fat
                dff.loc[:len(distance_art)-1, 'Distance to Arteriole'] = distance_art
                dff.loc[:len(distance_sinu)-1, 'Distance to Sinusoid'] = distance_sinu
                dff.to_csv('{}/{}.csv'.format(save_folder,sample_id), index=False)
                
                print('{}: {} out of {} : {} has been done'.format(comptype, i+1,len(zarr_),sample_id))
                
                del cluster_type
                del cell_centroid_
                del distance_bone
                del distance_fat
                del distance_art
                del distance_sinu
                
            else:
                raise ValueError("Only accepts [struct, celltype, cellneighbor].")
                
        elif comptype.split['_'][0] == 'celltype':
            
            cell_list = ['Adipo-MSC','B_cell', 'DC', 'Endothelial', 'Erythroid', 'GMP','Granulocyte/mast', 'HSPC',
                         'Macrophage', 'Megakaryocyte', 'Monocyte', 'Myeloid', 'Osteo-MSC', 'Plasma_cell', 
                         'SMC', 'Stromal', 'T_cell']
            
            cell_centroid = obj["new_cell_boundaries"]['geometry'].centroid
            cell_ids = adata.obs["cell_id_x"]
            cell_types = adata.obs["obj.anno_3_w_megs_w_stromal_x"]
            cell_bound = obj["new_cell_boundaries"]
            
            set_cell_ids = set(cell_ids)
            
            in_mask = [name in set_cell_ids for name in np.array(cell_centroid.index)]
            cell_centroid_ = cell_centroid[in_mask]
            
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
                
                x_data = np.zeros((len(cell_ids),4))
                x_data = x_data.astype(object)
                cell_types_temp = np.array(cell_types)
                x_data[:,0] = cell_types_temp
                x_data[:,1] = cell_centroid_.x
                x_data[:,2] = cell_centroid_.y
                x_data[:,3] = distance_cell
                header = ['Cell Type', 'x', 'y', 'Distance to Cell']
                dff = pd.DataFrame(x_data)
                dff.columns = header
                dff['Distance to Cell'] = dff['Distance to Cell'].astype(float)
                dff.set_index('Cell Type', inplace=True)
                dff.to_excel('{}/{}.xlsx'.format(subsavefolder,sample_id))
                print('{}: {} out of {} : {} --> {} has been done'.format(comptype, i+1,len(zarr_),sample_id, cell_to_name))
                del dff
                del distanec_cell
            
        elif comptype.split['_'][0] == 'cellneighbor':
            cluster_df = pd.read_csv('anndata_cc_3_n10.csv', index_col=0)

            cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']
            
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
                
                print('{}: {} out of {} : {} --> {} has been done'.format(comptype, i+1,len(zarr_),sample_id, cc))
                
                del dff
                del cluster_type_temp
        else:
            raise ValueError("Only accepts ['struct_struct', 'struct_celltype','struct_cellneighbor', 'celltype_celltype', 'cellneighbor_cellneighbor'].")
