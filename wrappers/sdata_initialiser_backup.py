import dask
dask.config.set({'dataframe.query-planning': False})


import sopa
#import spatialdata_xenium_explorer
from sopa.utils.data import uniform
import sopa.segmentation
import sopa.io
import spatialdata_io 
from spatialdata_io import xenium
import spatialdata as sd
import spatialdata_plot
import os
import matplotlib.pyplot as plt
import json
from shapely import MultiPolygon, Polygon, Point
from spatialdata.models import ShapesModel, TableModel
import geopandas as gpd
from spatialdata.transformations import (
    Affine,
    Identity,
    MapAxis,
    Scale,
    Sequence,
    Translation,
    get_transformation,
    get_transformation_between_coordinate_systems,
    set_transformation,
)
from spatialdata import transform
import pandas as pd

from spatialdata import SpatialData
from spatialdata.models import SpatialElement

import numpy as np
from sopa import segmentation
from spatialdata import polygon_query
import json

from geopandas import GeoDataFrame
from shapely.geometry import Polygon
import shapely

from shapely import LineString, Point
from shapely import get_coordinates

import anndata
import scanpy as sc

from spatialdata.models import TableModel
import matplotlib.pyplot as plt

import numpy as np
from spatialdata.models import ShapesModel

import math

import json

import geopandas

import anndata as ad


def get_intrinsic_cs(
    sdata: SpatialData, element: SpatialElement | str, name: str | None = None
) -> str:
    """Gets the name of the intrinsic coordinate system of an element

    Args:
        sdata: A SpatialData object
        element: `SpatialElement`, or its key
        name: Name to provide to the intrinsic coordinate system if not existing. By default, uses the element id.

    Returns:
        Name of the intrinsic coordinate system
    """
    if name is None:
        name = f"_{element if isinstance(element, str) else id(element)}_intrinsic"

    if isinstance(element, str):
        element = sdata[element]

    for cs, transform in get_transformation(element, get_all=True).items():
        if isinstance(transform, Identity):
            return cs

    set_transformation(element, Identity(), name)
    return name



def to_intrinsic(
    sdata: SpatialData, element: SpatialElement | str, element_cs: SpatialElement | str
) -> SpatialElement:
    """Transforms a `SpatialElement` into the intrinsic coordinate system of another `SpatialElement`

    Args:
        sdata: A SpatialData object
        element: `SpatialElement` to transform, or its key
        element_cs: `SpatialElement` of the target coordinate system, or its key

    Returns:
        The `SpatialElement` after transformation in the target coordinate system
    """
    if isinstance(element, str):
        element = sdata[element]
    cs = get_intrinsic_cs(sdata, element_cs)
    return sdata.transform_element_to_coordinate_system(element, cs)


class SDATA_INITIALISER:
    def __init__(self, keys_csv_path="/well/rittscher/users/qdv200/MPN/xenium/data/xenium_keys.csv", 
                 num=0, save_folder="/well/rittscher/users/qdv200/MPN/xenium/data/test_zarr_files/",
                 annotation_path = "/well/rittscher/users/qdv200/MPN/xenium/data/annotations/Xenium_annotations_DAN/",
                 cif_path = "/well/rittscher/users/qdv200/MPN/xenium/data/CIF/Coordinates_Scores/",
                 neg_annotation_save = "/well/rittscher/users/qdv200/MPN/xenium/data/transformed_negative_annotations/"
                 ):
        
        self.keys_csv = pd.read_csv(keys_csv_path)
        self.num = num
        #define which exp you want by xen_key
        self.sample_key = self.keys_csv["sample_key"][self.num]
        #get index of sample key
        self.sample_index = self.keys_csv[self.keys_csv["sample_key"] == self.sample_key].index[0]
        
        self.save_folder = save_folder

        self.base_path = self.keys_csv["base_path"][self.sample_index]
        self.exp_name = self.keys_csv["xen_exp"][self.sample_index]
        self.he_key = self.keys_csv["xenium file"][self.sample_index]
        self.he_image_path = "/well/rittscher/users/qdv200/MPN/xenium/data/morphology_files/ome_tif/" + self.he_key + ".ome.tif"


        self.meg_path = self.keys_csv["meg_path"][self.sample_index]
        self.meg_pheno_path = self.keys_csv["meg_pheno_path"][self.sample_index]

        self.annotation_path = annotation_path
        self.cif_path = cif_path

        self.mut_status = self.keys_csv["Mutation_status"][self.sample_index]
        self.run = self.keys_csv["Run"][self.sample_index]
        self.study_id = self.keys_csv["Study_ID"][self.sample_index]
        self.diagnosis = self.keys_csv["Diagnosis"][self.sample_index]
        self.diagnosis2 = self.keys_csv["Diagnosis_2"][self.sample_index]
        self.broad_diagnosis = self.keys_csv["Broad_diagnosis"][self.sample_index]

        self.neg_annotation_save = neg_annotation_save

        self.tissue_poly = self.keys_csv["tissue_poly"][self.sample_index]
        
    def initialise_sdata(self, save=True):
        sdata = xenium(os.path.join(self.base_path, self.exp_name))
        del sdata.images["morphology_mip"]
        
        if save == True:
            sdata.write(os.path.join(self.save_folder, self.sample_key + ".zarr"), overwrite=True)
        
        return sdata
    
    def initialise_sdata_from_zarr(self, zarr_path):
        
        sdata = sd.read_zarr(zarr_path)
        
        return sdata

class SDATA_CUSTOMISER(SDATA_INITIALISER):
    def __init__(self, keys_csv_path="/well/rittscher/users/qdv200/MPN/xenium/data/xenium_keys.csv", 
                 num=0, save_folder="/well/rittscher/users/qdv200/MPN/xenium/data/test_zarr_files/",
                 annotation_path = "/well/rittscher/users/qdv200/MPN/xenium/data/annotations/Xenium_annotations_DAN/",
                 cif_path = "/well/rittscher/users/qdv200/MPN/xenium/data/CIF/Coordinates_Scores/",
                 neg_annotation_save = "/well/rittscher/users/qdv200/MPN/xenium/data/transformed_negative_annotations/"
                 ):
        
        super().__init__(keys_csv_path, num, save_folder, annotation_path, cif_path, neg_annotation_save)
    
    
    def add_he_image(self, file_name=None):
        l_he_key = self.he_key.replace("x", "X")
        
        if file_name == None:
            zarr_path = os.path.join(self.save_folder, self.sample_key + '.zarr')
        else:
            zarr_path = file_name
        
        if print == True:
            print(f"sopa explorer add-aligned {zarr_path} {self.he_image_path} /well/rittscher/users/qdv200/MPN/xenium/data/landmarks_matrix/{l_he_key}_{self.sample_key}_matrix.csv --original-image-key morphology_focus")
        else:
            bashcommand = f"sopa explorer add-aligned {zarr_path} {self.he_image_path} /well/rittscher/users/qdv200/MPN/xenium/data/landmarks_matrix/{l_he_key}_{self.sample_key}_matrix.csv --original-image-key morphology_focus"
            return bashcommand
    
    #add in megakaryocytes
    
    def add_megakaryocytes(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        
        #get the megs out of the geojson file

        with open(self.meg_path) as f:
            data = json.load(f)

        megs = []
        for feature in data['features']:
            megs.append(feature['geometry']['coordinates'][0])

        megs_poly = [Polygon(meg) for meg in megs]
        
        #add phenotype to each megakaryocyte

        #load the phenotypes with pd

        meg_pheno = pd.read_csv(self.meg_pheno_path)

        #get the center point of each bounding box

        #split name at fromx

        megs_name = [meg.split("fromx")[1] for meg in meg_pheno['name']]

        meg_centroids = []

        for i in megs_name:
            xmin = int(i.split("_")[1])
            ymin = int(i.split("_")[3])
            xmax = int(i.split("_")[5])
            ymax = int(i.split("_")[7].replace(".png", ""))
            w = xmax - xmin
            h = ymax - ymin
            x = xmin + w/2
            y = ymin + h/2
            meg_centroids.append((x, y))

        #add the centroids back into df

        meg_pheno['centroid'] = meg_centroids
        
        gdf = gpd.GeoDataFrame({"geometry": megs_poly})
        
        meg_points = [Point(i) for i in meg_centroids]
        points_gdf = gpd.GeoDataFrame({"geometry": meg_points, "phenotype": meg_pheno['membership']})
        
        points_within = gpd.sjoin(gdf, points_gdf, how="left", predicate='contains')
        
        #put phenotype into gdf

        points_within.drop(columns=["index_right"], inplace=True)
        
        #remove duplicated indexes
        points_within.drop_duplicates(subset="geometry", inplace=True)
        
        new_gdf = ShapesModel.parse(points_within, transformations={"global": Identity()})
        
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
            
        transformed_gdf = transform(new_gdf, affine_he_to_xen, "global")
        
        sdata.shapes["transformed_megs"] = transformed_gdf
        
        set_transformation(sdata.shapes["transformed_megs"], Identity(), to_coordinate_system="global")
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_megs.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_megs.zarr"), overwrite=True)
        else:
        
            return sdata
        
    
    
    #change segmentation
    
    def add_megs_to_segmentation(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        print("Validating megs")
        
        meg_gdf = sdata.shapes["transformed_megs"]
        
        invalid_megs = meg_gdf.loc[~meg_gdf.geometry.is_valid]
        valid_megs = invalid_megs.make_valid()
        for i in list(valid_megs.index):
            #choose the biggest polygon from each 
            largest_poly = max(list(valid_megs[i].geoms), key=lambda a: a.area)
            invalid_megs.geometry[i] = largest_poly
        multi_megs = invalid_megs.loc[~invalid_megs.geometry.is_valid]
        print("number of multi megs", len(multi_megs))
        
        #now put the valid megs back into the sdata object

        for i in list(invalid_megs.index):
            meg_gdf.geometry[i] = invalid_megs.geometry[i]
        multi_megs = meg_gdf.loc[~meg_gdf.geometry.is_valid]
        
        #check only polygons in meg_gdf
        non_poly = meg_gdf[meg_gdf.geometry.type != "Polygon"]
        
        #take the largest polygon 

        for i in list(non_poly.index):
            print(i)
            #choose the biggest polygon from each 
            largest_poly = max(list(non_poly["geometry"][i].geoms), key=lambda a: a.area)
            #print(largest_poly)
            meg_gdf.geometry[i] = largest_poly
            
        #check for duplicates

        meg_gdf = meg_gdf.drop_duplicates(subset="geometry")
        
        sdata.shapes["transformed_megs"] = meg_gdf
        
        print("Validating and transforming cell boundaries")
        
        cell_gdf = sdata.shapes["cell_boundaries"]

        #check taht all the cells are valid

        invalid_cells = cell_gdf.loc[~cell_gdf.geometry.is_valid]
        
        valid_cells = invalid_cells.make_valid()
        
        #take the largest polygon 

        for i in list(valid_cells.index):
            #choose the biggest polygon from each 
            largest_poly = max(list(valid_cells[i].geoms), key=lambda a: a.area)
            invalid_cells.geometry[i] = largest_poly
            
        #check all rows in invalid_cells are polygons and not multipolygons

        non_poly = invalid_cells.loc[invalid_cells.geometry.geom_type != "Polygon"]
        
        #now put the invalid cells back into cell_gdf

        for i in list(invalid_cells.index):
            cell_gdf.geometry[i] = invalid_cells.geometry[i]
            
        #check all cells in cell_gdf are polygons and not multipolygons

        non_poly = cell_gdf.loc[cell_gdf.geometry.geom_type != "Polygon"]
        
        #take the largest polygon 

        for i in list(non_poly.index):
            print(i)
            #choose the biggest polygon from each 
            largest_poly = max(list(non_poly["geometry"][i].geoms), key=lambda a: a.area)
            #print(largest_poly)
            cell_gdf.geometry[i] = largest_poly
            
        #if I divide by microns are they in the same coordinate space

        # #divide each by 0.2125

        polys = cell_gdf.geometry.values

        #transform each poly to a list
        polys = [list(poly.exterior.coords) for poly in polys]

        #times each poly by 0.2125
        polys = [[(x/0.2125, y/0.2125) for x, y in poly] for poly in polys]
        polys
        
        from shapely.geometry import Polygon

        polys = [Polygon(poly) for poly in polys]
        cell_gdf["geometry"] = polys  
        
        set_transformation(cell_gdf, Identity(), "global")
        
        sdata.shapes["cell_boundaries"] = cell_gdf
        
        print("Add megs to cell boundaries")
        
        geo_df = sdata["transformed_megs"]
        old_geo_df = sdata["cell_boundaries"]

        #geo_df into same space as cells
        geo_df = to_intrinsic(sdata, geo_df, old_geo_df)  
        old_geo_df["cell_id"] = old_geo_df.index 
        
        gdf_join = old_geo_df.overlay(geo_df, how="difference")
        
        #check for multipolygons in gdf_join - these will be left over cells from the overlaps

        mp = gdf_join.loc[gdf_join.geometry.geom_type != "Polygon"]
        
        #explode each multipolygon into it's list of polygons
        mp = mp.explode()
        
        mp = mp.reset_index(drop=False)
        #make a new column by append the index to the cell_id
        mp["new_cell_id"] = mp["cell_id"] + "_" + mp["level_1"].astype(str)
        #make new_cell_id the index
        mp.index = list(mp["new_cell_id"])
        #drop some columns from mp
        mp = mp.drop(columns=["level_0", "level_1", "new_cell_id"])
        
        invalid = mp[~mp.geometry.is_valid]
        
        print("Check if polygons contain nuclei")
        
        nuclei = sdata.shapes["nucleus_boundaries"]
        
        # #divide each by 0.2125

        polys = nuclei.geometry.values

        #transform each poly to a list
        polys = [list(poly.exterior.coords) for poly in polys]

        #times each poly by 0.2125
        polys = [[(x/0.2125, y/0.2125) for x, y in poly] for poly in polys]
        
        #now convert back into a list of polygons

        from shapely.geometry import Polygon

        polys = [Polygon(poly) for poly in polys]

        #now convert back into a geodataframe
        nuclei["geometry"] = polys
        
        #now change teh transformation to identity because we have changed them manually

        set_transformation(nuclei, Identity(), "global")
        
        #add these back into sdata
        sdata["nucleus_boundaries"] = nuclei
        nuclei = sdata["nucleus_boundaries"]
        
        #get the centroid of each polygon
        centroids = [poly.centroid for poly in nuclei.geometry.values]

        #convert centroids into geodataframe

        df = pd.DataFrame(centroids, columns=["geometry"])
        df["radius"] = 1

        nuc_gdf = gpd.GeoDataFrame(
            df, geometry=df.geometry, 
        )
        
        #add index back in 
        nuc_gdf.index = list(nuclei.index)
        
        #remove multipolygons from new cell segmentation
        poly_gdf_join = gdf_join.loc[gdf_join.geometry.geom_type == "Polygon"]
        
        #get nuclei which are in megakaryocytes

        meg_gdf["nuclei"] = None

        for i in meg_gdf["geometry"].index:
            pol = meg_gdf["geometry"][i]
            #pol_gdf = gpd.GeoDataFrame({"geometry": [pol]})
            nuc_meg = nuc_gdf.within(pol)
            #get indexes where nuc_meg is true
            idx = list(nuc_meg[nuc_meg == True].index)
            meg_gdf["nuclei"][i] = idx
            
        #now remove the nuclei which are in megakaryocytes from nuc_gdf so that 

        nuc_in_megs = []

        for i in list(meg_gdf.index):
            nuc_in_megs.extend(meg_gdf["nuclei"][i])
            
        nuc_in_megs = list(set(nuc_in_megs))
        
        #remove these indexes from nuc_gdf

        sub_nuc_gdf = nuc_gdf.drop(nuc_in_megs)
        
        #so now if mp contains a nucleus put it back into cell_boundaries with a note to say altered

        mp["full_cell"] = False

        for i in mp["geometry"].index:
            pol = mp["geometry"][i]
            nuc_meg = nuc_gdf.within(pol)
            #get indexes where nuc_meg is true
            idx = list(nuc_meg[nuc_meg == True].index)
            if len(idx) > 0:
                mp["full_cell"][i] = True
        
        #get the indexes where they are the full cell

        full_cells = mp[mp["full_cell"] == True]
        
        #check taht all teh full cells are unique

        print("Full cells unique:", len(full_cells["cell_id"]) == len(list(set(full_cells["cell_id"]))))
        
        #drop full_cell column from full_cells

        full_cells = full_cells.drop(columns=["full_cell"])
        #add column
        full_cells["cell_status"] = "altered"
        
        #add full cells to poly_gdf_join

        poly_gdf_join["cell_status"] = "original"
        
        #concatenate full_cells and poly_gdf_join together

        full_cells_and_poly_gdf_join = pd.concat([poly_gdf_join, full_cells])
        
        print("Deal with small cells and artefacts")
        
        #now we need to subset mp to just the small cells and artefacts

        small_cells_and_artefacts = mp[mp["full_cell"] == False]
        
        #now I need to decide which are small cells and which are artefacts

        #get the size range of my cells in gdf_cells

        cell_area = cell_gdf.area
        
        art_area = small_cells_and_artefacts.area
        
        #define a small cell if less tahn half of median 
        #we will merge these into the cells
        small_cells = small_cells_and_artefacts[small_cells_and_artefacts.area < np.median(cell_area)/2]
        
        artefacts = small_cells_and_artefacts[small_cells_and_artefacts.area >= np.median(cell_area)/2]
        
        #find the nearest nucleus, that is not in a megakaryocyte, to each small cell
        #merge small cell into the cell that contains taht nucleus

        #not sure which point in polygon uses to calculate distance
        #relying on nuclei having same names as cells here

        nn_join = small_cells.sjoin_nearest(sub_nuc_gdf, max_distance=100, distance_col="distances")
        
        #okay now merge small cells into the cells in index_right

        cells_for_merge = list(set(nn_join["index_right"]))
        
        #check if each nuclei in cells_for_merge has a equivalent cell in index of full_cells_and_poly_gdf_join

        #check if list contains all elemnets of another list

        set(cells_for_merge).issubset(list(full_cells_and_poly_gdf_join["cell_id"]))
        
        print("Conglomerating all polygons together into new segmentation")
        
        small_cells_full_cells_and_poly_gdf_join = full_cells_and_poly_gdf_join.copy()

        for i in cells_for_merge:
            #get the polygons where index_right = i
            #get the small cells out 
            nn_sub = nn_join[nn_join["index_right"] == i]
            nn_sub.drop(columns=["distances", "radius", "index_right", "full_cell"], inplace=True)
            #print(nn_sub)
            
            #get the large cell taht it needs to be merged into
            nn_large = full_cells_and_poly_gdf_join[full_cells_and_poly_gdf_join["cell_id"] == i]
            #print(nn_large)
            
            #concatenate these two polygons into one geodataframe
            merged = pd.concat([nn_sub, nn_large])
            #print(merged)
            
            #now merge these polygons into one polygon using unary_union
            #expand so we can merge then de-expand - otherwise end up with a multipolygon
            union_poly = merged.buffer(10).unary_union.buffer(-10)
            #print(union_poly)
            
            small_cells_full_cells_and_poly_gdf_join.index = list(small_cells_full_cells_and_poly_gdf_join["cell_id"])
            small_cells_full_cells_and_poly_gdf_join["geometry"][i] = union_poly
            small_cells_full_cells_and_poly_gdf_join["cell_status"][i] = "altered_with_small_cells"
            
        print("add in artefacts and megakaryocytes")
        
        #okay now I just need to add in the artefact cells and megakaryocytes
        artefacts.drop(columns=["full_cell"], inplace=True)
        artefacts["cell_id"] = list(artefacts.index)
        artefacts["cell_status"] = "artefact"
        
        
        meg_gdf["cell_status"] = "megakaryocyte"
        meg_gdf.drop(columns="nuclei", inplace=True)
        meg_gdf["cell_id"] = list(meg_gdf.index)
        
        #now concatenate all together
        new_cell_gdf = pd.concat([small_cells_full_cells_and_poly_gdf_join, artefacts, meg_gdf])
        
        #if cell_status isn't original append new to cell_id and make index - else leave it so doesn't have to be updated

        for i in list(new_cell_gdf.index):
            if new_cell_gdf["cell_status"][i] != "original":
                new_cell_gdf["cell_id"][i] = str(i) + "_new"
                
        new_cell_gdf.index = list(new_cell_gdf["cell_id"])
        
        invalid = new_cell_gdf.loc[~new_cell_gdf.geometry.is_valid]
        
        #parse this with identity matrix
        set_transformation(new_cell_gdf, Identity(), "global")
        
        new_cell_gdf.index = new_cell_gdf["cell_id"]
    
        
        sdata.shapes["new_cell_boundaries"] = new_cell_gdf

        sdata["new_cell_boundaries"].assign(cell_status = list(new_cell_gdf["cell_status"]))
        
        #save this as a geodataframe
        
        metadata_save_folder = "/well/rittscher/users/qdv200/MPN/xenium/data/metadata_cell_status/"
        
        new_cell_gdf.to_csv(os.path.join(metadata_save_folder, self.sample_key + ".csv"), index=False)

        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_new_seg.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_new_seg.zarr"), overwrite=True)
        else:
        
            return sdata
        
        
    #aggregate with xenium_ranger
    
    def aggregate_with_xen_ranger(self, zarr_path=None, sdata=None, cell_poly_path="/well/rittscher/users/qdv200/MPN/xenium/data/new_cell_boundaries/"):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        cell_gdf = sdata["new_cell_boundaries"]
        
        #extract the largest polygon from each multipolygon

        for i in list(cell_gdf.index):
            #choose the biggest polygon from each 
            largest_poly = max(list(cell_gdf["geometry"][i].geoms), key=lambda a: a.area)
            #print(largest_poly)
            cell_gdf.geometry[i] = largest_poly
        
        
        #transform into a feature collection


        geojson = {
            "type": "FeatureCollection",
            "features": []
        }

        for i in list(cell_gdf.index):
            coords = list(cell_gdf["geometry"][i].exterior.coords)
            coords = [[x[0], x[1]] for x in coords]
            
            geojson["features"].append({
                "type": "Feature",
                "id": i,
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [coords]
                }
            })

        save_path = os.path.join(cell_poly_path, self.sample_key + ".geojson")

        with open(save_path, "w") as f:
            json.dump(geojson, f)
        
        output_folder = "/well/rittscher/users/qdv200/MPN/xenium/data/xen_ranger_outputs/"
        
        print(f"xeniumranger import-segmentation --id={os.path.join(output_folder, self.sample_key)} --xenium-bundle={os.path.join(self.base_path, self.exp_name)} --cells={save_path} --units=pixels --localcores=32 --localmem=128")
        

    #aggregate with sopa
    
    def aggregate_with_sopa(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        cell_gdf = sdata["new_cell_boundaries"]
        
        #extract the largest polygon from each multipolygon

        for i in list(cell_gdf.index):
            if cell_gdf["geometry"][i].geom_type == "MultiPolygon":
                #choose the biggest polygon from each 
                largest_poly = max(list(cell_gdf["geometry"][i].geoms), key=lambda a: a.area)
                #print(largest_poly)
                cell_gdf.geometry[i] = largest_poly
        
        sdata["new_cell_boundaries"] = cell_gdf
        
        new_cell_gdf = pd.read_csv(os.path.join("/well/rittscher/users/qdv200/MPN/xenium/data/metadata_cell_status/", self.sample_key + ".csv"))
        
        print(new_cell_gdf)
        
        old_anndata = sdata.tables["table"]
        
        #check for multipolygons in new_cell_boundaries

        non_poly = sdata["cell_boundaries"].loc[sdata["cell_boundaries"].geometry.geom_type != "Polygon"]
        print(len(non_poly))
        
        #https://gustaveroussy.github.io/sopa/tutorials/api_usage/
        aggregator = sopa.segmentation.Aggregator(sdata, image_key="morphology_focus", shapes_key="new_cell_boundaries", overwrite=False)
        
        #need to stop it changing the indices
        aggregator.compute_table(gene_column="feature_name", average_intensities=False)
        
        #change index (know this is fine because area the same)
        
        sdata["new_cell_boundaries"].index = list(new_cell_gdf["cell_id"])
        sdata.tables["table"].obs.index = list(new_cell_gdf["cell_id"])
        
        sdata.tables["table"].obs["cell_status"] = list(new_cell_gdf["cell_status"])
        
        sdata.tables["table"].obs["cell_id"] = list(new_cell_gdf["cell_id"])
        
        #subset annadat object so only genes

        sdata.tables["table"] = sdata.tables["table"][list(new_cell_gdf["cell_id"]), list(old_anndata.var_names)]
        
        # #now make region the cell boundaries so we can show cells to keep and remove
        # sdata.tables["table"].obs["region"] = "new_cell_boundaries"
        # sdata.set_table_annotates_spatialelement("table", region="new_cell_boundaries")
        
        categories = ["new_cell_boundaries"]
        n = len(sdata["table"])

        sdata["table"].obs["region"] = pd.Categorical(["new_cell_boundaries" for _ in range(n)], categories=categories)
        
        
        #now save
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_sopa_aggregated.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_sopa_aggregated.zarr"), overwrite=True)
        else:
            return sdata
        
    #add cif scores
    
    def add_512_cif_scores(self, zarr_path=None, sdata=None, save=False, cif_path="/well/rittscher/users/qdv200/MPN/xenium/data/CIF/Coordinates_Scores/", save_transformed_cif = False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        cif_score = pd.read_csv(os.path.join(cif_path, self.he_key + ".txt"), delim_whitespace = True)
        
        import shapely.wkt

        P = shapely.wkt.loads(self.tissue_poly)

        tiles_in_sample = []

        for i in cif_score.index:
            point = Point(cif_score["X"][i], cif_score["Y"][i])
            #print(point.within(P))
            if point.within(P) == True:
                tiles_in_sample.append(i)
        
        sub_cif_score = cif_score.iloc[tiles_in_sample]
        
        from shapely import geometry 
         
        tile_list = []  

        for i in list(sub_cif_score.index):
            minx = sub_cif_score["X"][i]
            miny = sub_cif_score["Y"][i]
            maxx = minx + 512
            maxy = miny + 512
            bbox_tuple = [minx, miny, maxx, maxy]
            bbox_polygon = geometry.box(*bbox_tuple)
            tile_list.append(bbox_polygon)
        
        cif_gdf = gpd.GeoDataFrame({"geometry": tile_list, "cif_score": sub_cif_score['Score']})
        
        cif_gdf = ShapesModel.parse(cif_gdf, transformations={"global": Identity()})
        
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        
        transformed_cif_gdf = transform(cif_gdf, affine_he_to_xen, "global")
        
        sdata.shapes["transformed_cif_tiles"] = transformed_cif_gdf
        
        set_transformation(sdata.shapes["transformed_cif_tiles"], Identity(), to_coordinate_system="global")
        
        if save_transformed_cif != False:
            transformed_cif_gdf.to_csv(os.path.join(save_transformed_cif, self.sample_key + "_512_transformed.csv"), index=False)
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_512_cif.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_512_cif.zarr"), overwrite=True)
        else:
            return sdata
        
    def add_1024_cif_scores(self, zarr_path=None, sdata=None, save=False, cif_path="/well/rittscher/users/qdv200/MPN/xenium/data/CIF/updated_CIF/", save_transformed_cif = False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        cif_score = pd.read_csv(os.path.join(cif_path, self.he_key + ".txt"), delim_whitespace = True)
        
        import shapely.wkt

        P = shapely.wkt.loads(self.tissue_poly)

        tiles_in_sample = []

        for i in cif_score.index:
            point = Point(cif_score["X"][i], cif_score["Y"][i])
            #print(point.within(P))
            if point.within(P) == True:
                tiles_in_sample.append(i)
        
        sub_cif_score = cif_score.iloc[tiles_in_sample]
        
        from shapely import geometry 
         
        tile_list = []  

        for i in list(sub_cif_score.index):
            minx = sub_cif_score["X"][i]
            miny = sub_cif_score["Y"][i]
            maxx = minx + 1024
            maxy = miny + 1024
            bbox_tuple = [minx, miny, maxx, maxy]
            bbox_polygon = geometry.box(*bbox_tuple)
            tile_list.append(bbox_polygon)
        
        cif_gdf = gpd.GeoDataFrame({"geometry": tile_list, "cif_score": sub_cif_score['Score']})
        
        cif_gdf = ShapesModel.parse(cif_gdf, transformations={"global": Identity()})
        
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        
        transformed_cif_gdf = transform(cif_gdf, affine_he_to_xen, "global")
        
        sdata.shapes["transformed_cif_1024_tiles"] = transformed_cif_gdf
        
        set_transformation(sdata.shapes["transformed_cif_1024_tiles"], Identity(), to_coordinate_system="global")
        
        if save_transformed_cif != False:
            transformed_cif_gdf.to_csv(os.path.join(save_transformed_cif, self.sample_key + "_1024_transformed.csv"), index=False)
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_1024_cif.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_1024_cif.zarr"), overwrite=True)
        else:
            return sdata
    
    
    def add_meg_phenotypes(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        #add in megakaryocytes from meg pipeline 

        #get the megs out of the geojson file

        with open(self.meg_path) as f:
            data = json.load(f)

        megs = []
        for feature in data['features']:
            megs.append(feature['geometry']['coordinates'][0])
            
        megs_poly = [Polygon(meg) for meg in megs]
        
        #add phenotype to each megakaryocyte

        #load the phenotypes with pd

        meg_pheno = pd.read_csv(self.meg_pheno_path)

        #get the center point of each bounding box

        #split name at fromx

        megs_name = [meg.split("fromx")[1] for meg in meg_pheno['name']]

        meg_centroids = []

        for i in megs_name:
            xmin = int(i.split("_")[1])
            ymin = int(i.split("_")[3])
            xmax = int(i.split("_")[5])
            ymax = int(i.split("_")[7].replace(".png", ""))
            w = xmax - xmin
            h = ymax - ymin
            x = xmin + w/2
            y = ymin + h/2
            meg_centroids.append((x, y))

        #add the centroids back into df

        meg_pheno['centroid'] = meg_centroids
        
        gdf = gpd.GeoDataFrame({"geometry": megs_poly})
        
        #make a points gdf

        meg_points = [Point(i) for i in meg_centroids]
        points_gdf = gpd.GeoDataFrame({"geometry": meg_points, "phenotype": meg_pheno['membership']})
        
        points_within = gpd.sjoin(gdf, points_gdf, how="left", predicate='contains')
        
        #put phenotype into gdf

        points_within.drop(columns=["index_right"], inplace=True)
        
        #remove duplicated indexes
        points_within.drop_duplicates(subset="geometry", inplace=True)
        
        new_gdf = ShapesModel.parse(points_within, transformations={"global": Identity()})
        
        sdata["transformed_megs"]["phenotype"] = list(new_gdf["phenotype"])
        
        #append new to each of the indexes

        meg_gdf = sdata["transformed_megs"]
        meg_gdf["cell_id"] = list(meg_gdf.index)

        for i in list(meg_gdf.index):
            meg_gdf["cell_id"][i] = str(i) + "_new"
            
        meg_gdf.index = list(meg_gdf["cell_id"])
        
        sdata.tables["table"].obs["meg_phenotype"] = "non_meg"

        for i in list(meg_gdf.index):
            sdata.tables["table"].obs["meg_phenotype"][i] = meg_gdf["phenotype"][i]
            
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_meg_pheno.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_meg_pheno.zarr"), overwrite=True)
        else:
            return sdata
    
    def add_negative_annotations(self, zarr_path=None, sdata=None, save=False, save_neg_annotations="/well/rittscher/users/qdv200/MPN/xenium/data/transformed_negative_annotations/"):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        #now add in the negative and positive annotations
        #visualise dan's annotations

        ann_path = os.path.join(self.annotation_path, self.sample_key)
        files = os.listdir(ann_path)
        
        #import the negative selections
        neg_polys = pd.read_csv(os.path.join(ann_path, self.sample_key + "_N_C.csv"), sep=",", comment='#')

        #make index a column
        neg_polys = neg_polys.reset_index()
        #make the first row the header
        neg_polys.columns = neg_polys.iloc[0]
        neg_polys = neg_polys[1:]
        neg_polys["X"] = neg_polys["X"].astype(float).div(0.2125)
        neg_polys["Y"] = neg_polys["Y"].astype(float).div(0.2125)

        neg_uni_selection = neg_polys["Selection"].unique()
        
        from shapely import MultiPolygon, Polygon

        crushed_polygons = []
        crushed_polygons_num = []

        for i in range(len(neg_uni_selection)):
            
            neg_s = neg_polys[neg_polys["Selection"] == neg_uni_selection[i]]
            poly = Polygon(neg_s[["X", "Y"]].values.tolist())
            crushed_polygons.append(poly)
            
            #get out the number and add it to a list, to make a column in geodataframe
            num = neg_uni_selection[i].split(" ")[-1]
            crushed_polygons_num.append(num)
        
        #add the intratrebecular regions into the sdata object

        gdf2 = gpd.GeoDataFrame({"geometry": crushed_polygons, "number": crushed_polygons_num})
        gdf2 = ShapesModel.parse(gdf2, transformations={"global": Identity()})
        sdata.shapes["crushed_regions"] = gdf2
        
        #get the he transformed version to save for Hosuk

        affine_glob_to_he = get_transformation_between_coordinate_systems(sdata, "global", sdata.images[self.he_key])
        affine_glob_to_he

        transformed_gdf = transform(gdf2, affine_glob_to_he, "he_space")
        
        if save_neg_annotations != False:
            transformed_gdf.to_csv(os.path.join(save_neg_annotations, self.sample_key + "_N.csv"), index=False)
        
        #now find the cells to remove under the negative polygons
        
        #create a faster way of doing this 

        sdata.tables["table"].obs["annotation"] = "keep"

        cells_to_remove = []

        for i in list(sdata["crushed_regions"].index):
            poly = sdata["crushed_regions"].geometry[i]
            outcome = sdata["new_cell_boundaries"].within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            cells_to_remove.extend(cell_ids)
        
        for i in cells_to_remove:
            sdata.tables["table"].obs["annotation"][i] = "remove"
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_neg_poly.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_neg_poly.zarr"), overwrite=True)
        else:
            return sdata
        
    def add_positive_annotations(self, zarr_path=None, sdata=None, save=False, save_pos_annotations="/well/rittscher/users/qdv200/MPN/xenium/data/transformed_positive_annotations/"):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata   
            
        #import the positive selections
        
        #now add in the negative and positive annotations
        #visualise dan's annotations

        ann_path = os.path.join(self.annotation_path, self.sample_key)
        files = os.listdir(ann_path)
        pos_polys = pd.read_csv(os.path.join(ann_path, self.sample_key + "_P_C.csv"), sep=",", comment='#')

        #make index a column
        pos_polys = pos_polys.reset_index()

        #make the first row the header
        pos_polys.columns = pos_polys.iloc[0]
        pos_polys = pos_polys[1:]

        #each pixel is 0.2125 microns - xenium workflow but can also see in image metadata
        #scale each column
        pos_polys["X"] = pos_polys["X"].astype(float).div(0.2125)
        pos_polys["Y"] = pos_polys["Y"].astype(float).div(0.2125)

        #get unique values of selection so can iterate through
        pos_uni_selection = pos_polys["Selection"].unique() 
        
        #add the positive polygons as a shape to the sdata so that they are saved

        from shapely import MultiPolygon, Polygon

        intra_polygons = []
        intra_polygons_num = []

        for i in range(len(pos_uni_selection)):
            
            pos_s = pos_polys[pos_polys["Selection"] == pos_uni_selection[i]]
            poly = Polygon(pos_s[["X", "Y"]].values.tolist())
            intra_polygons.append(poly)
            
            #get out the number and add it to a list, to make a column in geodataframe
            num = pos_uni_selection[i].split(" ")[-1]
            intra_polygons_num.append(num)
            
        #add the intratrebecular regions into the sdata object

        gdf2 = gpd.GeoDataFrame({"geometry": intra_polygons, "number": intra_polygons_num})
        gdf2 = ShapesModel.parse(gdf2, transformations={"global": Identity()})
        sdata.shapes["intertrabecular_regions"] = gdf2
        
        #get the he transformed version to save for Hosuk

        affine_glob_to_he = get_transformation_between_coordinate_systems(sdata, "global", sdata.images[self.he_key])
        affine_glob_to_he

        transformed_gdf = transform(gdf2, affine_glob_to_he, "he_space")
        
        if save != False:
            transformed_gdf.to_csv(os.path.join(save_pos_annotations, self.sample_key + "_P.csv"), index=False)

        #add in positive intratrebecular annotations

        #make a column to of cells to remove based on annotations

        sdata.tables["table"].obs["it_regions"] = "non_intertrabecular"

        it_cells = []
        cell_ids_it = []

        for i in list(sdata["intertrabecular_regions"].index):
            print(i)
            poly = sdata["intertrabecular_regions"].geometry[i]
            outcome = sdata["new_cell_boundaries"].within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            cell_ids_it.extend(cell_ids)
            it = [i]*len(cell_ids)
            it_cells.extend(it) 
        
        #make a dataframe

        df = pd.DataFrame(it_cells, columns=["it_regions"])
        df.index = cell_ids_it
        
        df2 = df[~df.index.duplicated(keep='first')]
        
        sdata.tables["table"].obs.update(df2)
        
        print(set(sdata.tables["table"].obs["it_regions"]))

        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_pos_poly.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_pos_poly.zarr"), overwrite=True)
        else:
            return sdata 
        
    
    def add_metadata(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        #add in the metadata
        
        #append sample name to indexes of new_cell_boundaries and cell_id and index of obs table

        sdata.tables["table"].obs["cell_id"] = sdata.tables["table"].obs["cell_id"].str.replace("_new", "")
        sdata.tables["table"].obs["cell_id"] = sdata.tables["table"].obs["cell_id"] + "_" + self.sample_key
        
        
        sdata["new_cell_boundaries"].index = list(sdata.tables["table"].obs["cell_id"])
        
        sdata.tables["table"].obs.index = list(sdata.tables["table"].obs["cell_id"])
        
        #now add in the other column 

        sdata.tables["table"].obs["mutation_status"] = self.mut_status
        sdata.tables["table"].obs["run"] = self.run
        sdata.tables["table"].obs["study_id"] = self.study_id
        sdata.tables["table"].obs["diagnosis"] = self.diagnosis
        sdata.tables["table"].obs["diagnosis2"] = self.diagnosis2
        sdata.tables["table"].obs["broad_diagnosis"] = self.broad_diagnosis
        
        sdata.tables["table"].obs["meg_phenotype"] = sdata.tables["table"].obs["meg_phenotype"].astype(str)
        
        sdata.tables["table"].obs["it_regions"] = sdata.tables["table"].obs["it_regions"].astype(str)
        
        
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_metadata.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_metadata.zarr"), overwrite=True)
        else:
            return sdata
        
    
    def add_cell_annotations(self, zarr_path=None, sdata=None, save=False, anno_path = "/well/rittscher/users/qdv200/MPN/xenium/data/object_metadata/slimmed_all_metadata_03.csv", original_cluster_name="SCT_snn_res.0.3"):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        #add in the cell annotations
        
        anno = pd.read_csv(anno_path)
        
        #make sample id column by extracting last 8 characters from cell_id

        anno["sample_id"] = anno["cell_id"].str[-8:]
        
        anno_sub = anno[anno["sample_id"] == self.sample_key]
        
        anno_sub = anno_sub[["cell_id", original_cluster_name, "clusters"]]
        
        sdata.tables["table"].obs = pd.merge(sdata.tables["table"].obs, anno_sub, left_on="cell_id", right_on="cell_id", how="left")
        
        #now change nan (remove so not been annottaed) to remove rather than nan so can save
        
        sdata.tables["table"].obs["clusters"] = sdata.tables["table"].obs["clusters"].where(sdata.tables["table"].obs["annotation"] != "remove", "remove")
        sdata.tables["table"].obs["clusters"] = sdata.tables["table"].obs["clusters"].where(sdata.tables["table"].obs["clusters"].isnull() != True , "remove")
        
        #now add megs in from morphology
        sdata.tables["table"].obs["clusters_w_megs"] = sdata.tables["table"].obs["clusters"].where(sdata.tables["table"].obs["cell_status"] != "megakaryocyte" , "Megs")
        
        sdata.tables["table"].obs.index = list(sdata.tables["table"].obs["cell_id"])
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_annotations.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_annotations.zarr"), overwrite=True)
        else:
            return sdata
        
    def change_cell_annotations(self, zarr_path=None, sdata=None, save=False, map_path = "/well/rittscher/users/qdv200/MPN/xenium/data/cell_annotation/revised_feature_map.csv", original_cluster_name="SCT_snn_res.0.3"):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        df_map = pd.read_csv(map_path)  
        map_dict = dict(zip(df_map[original_cluster_name], df_map["cell_type"]))
        sdata.tables["table"].obs["new_clusters"] =  sdata.tables["table"].obs[original_cluster_name]
        sdata.tables["table"].obs["new_clusters"] = sdata.tables["table"].obs["new_clusters"].map(map_dict)
        sdata.tables["table"].obs["new_clusters"] = sdata.tables["table"].obs["new_clusters"].where(sdata.tables["table"].obs["annotation"] != "remove", "remove")
        sdata.tables["table"].obs["new_clusters"] = sdata.tables["table"].obs["new_clusters"].where(sdata.tables["table"].obs["new_clusters"].isnull() != True , "remove")
                
        #now add megs in from morphology
        sdata.tables["table"].obs["new_clusters_w_megs"] = sdata.tables["table"].obs["new_clusters"].where(sdata.tables["table"].obs["cell_status"] != "megakaryocyte" , "Mks")
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_changed_annotations.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_changed_annotations.zarr"), overwrite=True)
        else:
            return sdata 
        
    def add_cif_abundance(self, zarr_path=None, sdata=None, save=False, tile_path = "/well/rittscher/users/qdv200/MPN/xenium/data/score_tiles/transformed_512_cif/", tile_size = 512, cell_polys="new_cell_boundaries_with_adipocytes"):

        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        
        df_cif = pd.read_csv(os.path.join(tile_path, self.sample_key + "_" + str(tile_size) + "_transformed.csv"))
        
        gdf = GeoDataFrame(df_cif, geometry=list(df_cif["geometry"].apply(shapely.wkt.loads)))
        
        if tile_size == 512:
            shape_name = "transformed_cif_tiles"
            circle_name = "transformed_cif_circles"
        else:
            shape_name = "transformed_cif_1024_tiles"
            circle_name = "transformed_cif_1024_circles"
        
        if all(sdata[shape_name].geom_equals(gdf, align=False)) == False:
            print("Geometry does not match")
            return
        else:
            print("Geometry matches")
        
        
        new_index = [str(i) + "_" + self.sample_key for i in range(len(gdf))]
        df_cif.index = new_index
        sdata[shape_name].index = new_index
        
        centers = gdf["geometry"].centroid
        df_cif["center"] = list(centers)
        
        x_center = [get_coordinates(i).tolist()[0][0] for i in df_cif["center"]]
        y_center = [get_coordinates(i).tolist()[0][1] for i in df_cif["center"]]
        
        df_cif["center_x"] = x_center
        df_cif["center_y"] = y_center
        df_cif.drop(columns=["center"], inplace=True)
        df_cif.drop(columns=["geometry"], inplace=True)
        
        cats = list(sdata.tables["table"].obs["new_clusters_w_megs"].unique())
        
        df_X = pd.DataFrame(index=new_index, columns=cats)
        
        for cat in cats:
            df_X[cat + "_prop"] = None

        df_X["ncells"] = None
            

        for i in list(sdata[shape_name].index):
            #print(i)
            poly = sdata[shape_name].geometry[i]
            outcome = sdata["new_cell_boundaries"].within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            #print(cell_ids)
            #now extract the new_clusters_w_megs associated with these cells
            counts = sdata.tables["table"].obs.loc[cell_ids, "new_clusters_w_megs"].value_counts()
            counts = counts.to_dict()
            
            df_X.loc[i] = counts
            
            
            #get sum of all values in dict
            total = sum(counts.values())
            
            #add keys to counts
        
            for cat in cats:
                if cat in counts.keys():
                    if total == 0:
                        df_X.loc[i, cat + "_prop"] = 0
                    else:
                        df_X.loc[i, cat + "_prop"] = counts[cat] / total
                else:
                    df_X.loc[i, cat + "_prop"] = None

            df_X.loc[i, "ncells"] = total 
        
        df_X.fillna(0, inplace=True)
        
        df_cif["annotation"] = "keep"

        for i in df_X.index:
            if df_X.loc[i, "ncells"] == 0 or df_X.loc[i, "remove_prop"] >= 0.3:
                df_cif.loc[i, "annotation"] = "remove"
                
        df_cif["region"] = shape_name
        df_cif["instance_id"] = list(df_cif.index)
        
        adata = anndata.AnnData(X=df_X, obs=df_cif)
        
        table = TableModel.parse(adata, region=shape_name, region_key="region", instance_key="instance_id")
        
        
        
        table_name = "cif_" + str(tile_size) + "_table"
        
        sdata.tables[table_name] = table
        
        arr = np.array(list(sdata.tables[table_name].X), dtype=float)
        sdata.tables[table_name].X = arr
        
        #now add cif circles
        
        #get the size of the cif tiles

        x, y = sdata[shape_name]["geometry"][0].exterior.coords.xy

        x = x.tolist()
        y = y.tolist()

        max_radii = (x[1] - x[0])/2
        print(max_radii)
        
        gdf = sdata[shape_name].copy()
        gdf["geometry"] = gdf["geometry"].centroid
        
        gdf["radius"] = 0
        for i in gdf.index:
            cs = sdata.tables[table_name].obs["cif_score"][i]
            
            #if cs is more than 0 and isn't nan
            if cs > 0 or math.isnan(cs) == False:
                new_radii = max_radii * cs
            else:
                new_radii = np.nan
            
            if new_radii == 0:
                new_radii = 0.1
                
            gdf.loc[i, "radius"] = new_radii
        
        gdf["radius"].fillna(0.1, inplace=True)
        
        print(min(gdf["radius"]))
        
        new_gdf = ShapesModel.parse(gdf)
        
        #print(new_gdf)
        
        sdata[circle_name] = new_gdf
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_" + str(tile_size) + "_cif_abundance.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_" + str(tile_size) + "_cif_abundance.zarr"), overwrite=True)
        else:
            return sdata 
        
        
    def add_poly_abundance(self, zarr_path=None, sdata=None, save=False):

        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        shape_name = "intertrabecular_regions"

        df = sdata.tables["table"].obs

        gdf = sdata["intertrabecular_regions"]

        old_index = list(gdf.index)
        old_index = map(str, old_index)

        new_index = [self.sample_key + "_R" + str(i) for i in old_index]
        gdf.index = new_index
        sdata["intertrabecular_regions"].index = new_index

        df_cif = pd.DataFrame(index=list(gdf.index), columns=["polygon_id"])


        centers = gdf["geometry"].centroid
        df_cif["center"] = list(centers)

        x_center = [get_coordinates(i).tolist()[0][0] for i in df_cif["center"]]
        y_center = [get_coordinates(i).tolist()[0][1] for i in df_cif["center"]]

        df_cif["center_x"] = x_center
        df_cif["center_y"] = y_center
        df_cif.drop(columns=["center"], inplace=True)
        df_cif["polygon_id"] = list(df_cif.index)

        cats = list(sdata.tables["table"].obs["new_clusters_w_megs"].unique())


        sub = sdata.tables["table"].obs

        #subset to only index and it_regions

        sub = sub[["it_regions", "cell_id"]]

        sub = sub[sub["it_regions"] != "non_intertrabecular"]

        #change sub 

        #sub["sample_id"] = sub["cell_id"].str[-8:]

        sub.drop(columns=["cell_id"], inplace=True)

        concat_func = lambda x,y: x + "_R" + y

        it_regions = list(map(float, sub["it_regions"]))
        it_regions = list(map(int, it_regions))
        it_regions = list(map(str, it_regions))

        sub["it_regions"] = it_regions

        region_name = list(map(concat_func, [self.sample_key]*len(sub["it_regions"]), sub["it_regions"]))

        #merge back into the main dataframe

        sub["region_name"] = region_name

        sub = sub[["region_name"]]

        df_X = df.merge(sub, how="left", left_index=True, right_index=True)

        df_X["region_name"].unique()

        df_X["region_name"].fillna("non_intertrabecular", inplace=True)

        sdata.tables["table"].obs["region_name"] = df_X["region_name"]

        region_adata = sdata.tables["table"]

        cats = list(region_adata.obs["new_clusters_w_megs"].unique())
        sample_id = list(region_adata.obs["region_name"].unique())

        sample_id = [i for i in sample_id if i != "non_intertrabecular"]
        
        df_X = pd.DataFrame(index=list(df_cif.index), columns=cats)

        for cat in cats:
            df_X[cat + "_prop"] = None

        df_X["ncells"] = None
        
        df_cif["region"] = shape_name

        df_cif["instance_id"] = df_cif.index


        for i in list(df_X.index):
            print(i)

            sub = region_adata.obs[region_adata.obs["region_name"] == i]
            cell_ids = list(sub.index)

            counts = region_adata.obs.loc[cell_ids, "new_clusters_w_megs"].value_counts()
            counts = counts.to_dict()
            
            df_X.loc[i] = counts
            
            #get sum of all values in dict
            total = sum(counts.values())
            
            #add keys to counts

            for cat in cats:
                if cat in counts.keys():
                    if total == 0:
                        df_X.loc[i, cat + "_prop"] = 0
                    else:
                        df_X.loc[i, cat + "_prop"] = counts[cat] / total
                else:
                    df_X.loc[i, cat + "_prop"] = None

            df_X.loc[i, "ncells"] = total 

        df_X.fillna(0, inplace=True)

        df_cif["diagnosis2"] = [self.diagnosis] * len(df_cif.index)
            
        df_cif["sample_id"] = df_cif.index.str[0:8]

        adata = anndata.AnnData(X=df_X, obs=df_cif)
                
        table = TableModel.parse(adata, region=shape_name, region_key="region", instance_key="instance_id")

        table_name = shape_name + "_table"

        sdata.tables[table_name] = table


        arr = np.array(list(sdata.tables[table_name].X), dtype=float)
        sdata.tables[table_name].X = arr

        #now add cif circles

        #get the size of the cif tiles

        gdf = sdata[shape_name].copy()
        gdf["geometry"] = gdf["geometry"].centroid

        gdf["radius"] = 3000

        new_gdf = ShapesModel.parse(gdf)

        #print(new_gdf)

        circle_name = shape_name + "_circles"

        sdata[circle_name] = new_gdf
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_poly_abundance.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_poly_abundance.zarr"), overwrite=True)
        else:
            return sdata 
        
    def add_bone_segmentation(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
        
        bone_seg_path = "/well/rittscher/users/qdv200/MPN/xenium/data/bone_segmentation/"

        bone_geo_path = bone_seg_path + self.sample_key + ".geojson"

        with open(bone_geo_path) as f:
            gj = json.load(f)
            
        gdf = GeoDataFrame.from_features(gj)
        
        new_gdf = ShapesModel.parse(gdf, transformations={"global": Identity()})
        
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        
        transformed_gdf = transform(new_gdf, affine_he_to_xen, "global")
        
        transformed_gdf.index = [f"{self.sample_key}_B{i}" for i in range(len(transformed_gdf))]
        
        sdata.shapes["transformed_bone"] = transformed_gdf
        
        set_transformation(sdata.shapes["transformed_bone"], Identity(), to_coordinate_system="global")
            
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_bone.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_bone.zarr"), overwrite=True)
        else:
            return sdata
    
    
    def add_adipocytes(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
    
        base = "/well/rittscher/users/qdv200/MPN/xenium/data/adipocytes"
        
        with open(os.path.join(base, self.sample_key + ".geojson")) as f:
            data = json.load(f)
        
        geodf = GeoDataFrame.from_features(data["features"], columns=["geometry"])
        
        geodf = geodf[geodf["geometry"].geom_type != "LineString"]
            
        geodf.reset_index(drop=True, inplace=True)
        
        new_gdf = ShapesModel.parse(geodf, transformations={"global": Identity()})
        
        sdata.shapes["adipocytes"] = new_gdf
        
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        
        transformed_gdf = transform(new_gdf, affine_he_to_xen, "global")
        
        new_transformed_gdf = transformed_gdf.copy()
        new_transformed_representation = transformed_gdf.copy()

        #new_transformed_gdf["geometry"] = transformed_gdf["geometry"].apply(lambda x: x.buffer(-50) if x.area > 8000 else x.buffer(-30))
        
        centers = new_transformed_representation["geometry"].centroid
        #now make each of those centers into a circle of radius 30
        new_transformed_representation["geometry"] = [Point(i) for i in centers]
        new_transformed_representation["geometry"] = new_transformed_representation["geometry"].apply(lambda x: x.buffer(15))
        
        new_transformed_representation.index = [f"{self.sample_key}_A{i}" for i in range(len(new_transformed_representation))]

        new_transformed_gdf.index = [f"{self.sample_key}_A{i}" for i in range(len(new_transformed_gdf))]
        
        sdata.shapes["transformed_adipocytes"] = new_transformed_gdf
        
        sdata.shapes["transformed_adipocytes_circles"] = new_transformed_representation
            
        set_transformation(sdata.shapes["transformed_adipocytes"], Identity(), to_coordinate_system="global")
        set_transformation(sdata.shapes["transformed_adipocytes_circles"], Identity(), to_coordinate_system="global")
        
        df1 = sdata["new_cell_boundaries"]
        df2 = sdata["transformed_adipocytes_circles"]

        res_symdiff = geopandas.overlay(df1, df2, how='difference')
        
        sdata["new_cell_boundaries_without_adipocytes"] = res_symdiff
        
        sdata["new_cell_boundaries_without_adipocytes"].index = [f"{self.sample_key}_{i}" for i in range(len(sdata["new_cell_boundaries_without_adipocytes"]))]
        
        sdata["new_cell_boundaries_with_adipocytes"] = pd.concat([sdata["new_cell_boundaries_without_adipocytes"], sdata["transformed_adipocytes_circles"]])

        sdata["new_cell_boundaries_with_full_adipocytes"] = pd.concat([sdata["new_cell_boundaries_without_adipocytes"], sdata["transformed_adipocytes"]])
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_unaggregated_adipocytes.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_unaggregated_adipocytes.zarr"), overwrite=True)
        else:
            return sdata
    
    def add_adipo_table(self, zarr_path=None, sdata=None, save=False):
    
    
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        #change index
        gdf = sdata["new_cell_boundaries_with_adipocytes"]
        fat_indexes = gdf.index[gdf.index.str.contains("_A")].tolist()
        gdf_no_fat = gdf[~gdf.index.isin(fat_indexes)]
        cell_gdf = sdata["new_cell_boundaries"]
        cell_indexes = list(cell_gdf.index)
        all_indexes = cell_indexes + fat_indexes
        sdata["new_cell_boundaries_with_adipocytes"].index = all_indexes
        sdata["new_cell_boundaries_with_full_adipocytes"].index = all_indexes
        
        #make anndata of fat
        
        adata = sdata.tables["table"]
        col_names = adata.obs.columns
        fat_df = pd.DataFrame(index=fat_indexes, columns=col_names)
        fat_df["region"] = "new_cell_boundaries"
        fat_df["slide"] = "morphology_focus"
        fat_df["cell_id"] = fat_indexes
        fat_df["area"] = None
        fat_df["cell_status"] = "adipocyte"
        fat_df["meg_phenotype"] = "non_meg"
        fat_df["annotation"] = None
        fat_df["it_regions"] = None
        fat_df["mutation_status"] = self.mut_status
        fat_df["run"] = self.run
        fat_df["study_id"] = self.study_id
        fat_df["diagnosis"] = self.diagnosis
        fat_df["diagnosis2"] = self.diagnosis2
        fat_df["broad_diagnosis"] = self.broad_diagnosis
        fat_df["SCT_snn_res.0.3"] = "adipocyte"
        fat_df["clusters"] = "adipocyte"
        fat_df["clusters_w_megs"] = "adipocyte"
        fat_df["new_clusters_w_megs"] = "adipocyte"
        fat_df["new_clusters"] = "adipocyte"
        fat_df["region_name"] = None
        
        #now create the X matrix
        from scipy.sparse import csr_matrix
        counts = csr_matrix(np.zeros(shape=(len(fat_indexes), len(adata.var_names))), dtype=np.float32)
        fat_adata = ad.AnnData(counts)
        
        fat_adata.obs_names = fat_indexes
        fat_adata.var_names = list(adata.var_names)
        
        fat_adata.obs = fat_df
        
        
        new_adata = ad.concat([adata, fat_adata], join="outer")
        
        new_adata.uns = adata.uns
        
        sdata_table = TableModel.parse(new_adata)
        
        sdata.tables["table"] = sdata_table
        
        sdata.tables["table"].obs["region"] = "new_cell_boundaries_with_adipocytes"
        sdata.set_table_annotates_spatialelement("table", region="new_cell_boundaries_with_adipocytes")
        
        sdata.tables["table"].obs["region"] = sdata.tables["table"].obs["region"].astype("category")
        sdata.tables["table"].obs["it_regions"] = sdata.tables["table"].obs["it_regions"].astype("str")
        sdata.tables["table"].obs["it_regions"] = sdata.tables["table"].obs["it_regions"].astype("category")
        sdata.tables["table"].obs["SCT_snn_res.0.3"] = sdata.tables["table"].obs["SCT_snn_res.0.3"].astype("str")
        sdata.tables["table"].obs["SCT_snn_res.0.3"] = sdata.tables["table"].obs["SCT_snn_res.0.3"].astype("category")
                
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_adipo_table.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_adipo_table.zarr"), overwrite=True)
        else:
            return sdata
    
    def distance_from_bone(self, zarr_path=None, sdata=None, save=False):
    
    
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        sdata["new_cell_boundaries_with_full_adipocytes"].index = sdata["new_cell_boundaries_with_adipocytes"].index.tolist()
        
        bones = sdata["transformed_bone"]
        cells = sdata["new_cell_boundaries"]

        cells_w_bones = geopandas.sjoin_nearest(cells, bones, distance_col="distances", how="left", max_distance=20000)
        
        def diff(first, second):
            second = set(second)
            return [item for item in first if item not in second]

        #get duplicates in a list python

        def get_duplicates(lst):
            seen = set()
            seen_add = seen.add
            # adds all elements it doesn't know yet to seen and all other to seen_twice
            seen_twice = set( x for x in lst if x in seen or seen_add(x) )
            # turn the set into a list (as requested)
            return list( seen_twice )

        #get the indexes of the duplicates

        duplicates = get_duplicates(cells_w_bones.index)
        
        df_dropped = cells_w_bones.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        df_dropped = df_dropped.rename(columns={"index_right": "nearest_bone", "distances": "distance_to_bone"})

        df_dropped = df_dropped[["nearest_bone", "distance_to_bone"]]
        df = sdata.tables["table"].obs.merge(df_dropped, left_index=True, right_index=True)

        sdata.tables["table"].obs["nearest_bone"] = df["nearest_bone"]
        sdata.tables["table"].obs["distance_to_bone"] = df["distance_to_bone"]
        
        sdata.tables["table"].obs["bone_region"] = "intertrabecular"

        sdata.tables["table"].obs.loc[(sdata.tables["table"].obs["distance_to_bone"] > 0) & (sdata.tables["table"].obs["distance_to_bone"] < 300), "bone_region"] = "peritrabecular"

        #find the cells which are inside a bone polygon

        within_bone = []

        for i in list(bones.index):
            poly = bones.geometry[i]
            outcome = sdata["new_cell_boundaries"].within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            within_bone.extend(cell_ids)

        #for each cell where distance to bone is 0 call it endosteal

        sdata.tables["table"].obs.loc[sdata.tables["table"].obs["distance_to_bone"] < 1, "bone_region"] = "endosteal"

        sdata.tables["table"].obs.loc[within_bone, "bone_region"] = "within_bone"
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_bone_distance.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_bone_distance.zarr"), overwrite=True)
        else:
            return sdata
        
    
    def fibrotic_distance(self, zarr_path=None, sdata=None, save=False):
        #categorise each tile into a category
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata

        sdata.tables["cif_512_table"].obs["cif_category"] = "Unknown"

        sdata.tables["cif_512_table"].obs.loc[sdata.tables["cif_512_table"].obs["cif_score"] >= 0, "cif_category"] = "MF-0"

        sdata.tables["cif_512_table"].obs.loc[sdata.tables["cif_512_table"].obs["cif_score"] > 0.25, "cif_category"] = "MF-1"

        sdata.tables["cif_512_table"].obs.loc[sdata.tables["cif_512_table"].obs["cif_score"] > 0.5, "cif_category"] = "MF-2"

        sdata.tables["cif_512_table"].obs.loc[sdata.tables["cif_512_table"].obs["cif_score"] > 0.75, "cif_category"] = "MF-3"
        
        df = sdata.tables["cif_512_table"].obs[sdata.tables["cif_512_table"].obs["cif_category"] == "MF-3"]
        
        high_cif = list(df.index)
        gdf = sdata["transformed_cif_tiles"].loc[high_cif]
        gdf["geometry"] = gdf.buffer(1)
        combined_polygons = gdf.dissolve().explode()
        #make a new index unique for each niche and remove starting 0 column index
        combined_polygons = combined_polygons.reset_index(drop=True)
        #append sample_key and MF_niche to the index
        combined_polygons["index"] = self.sample_key + "_MF_niche_" + combined_polygons.index.astype(str)
        #make index the index
        combined_polygons.index = list(combined_polygons["index"])

        #get the index that of MF that are in each niche

        sdata.tables["cif_512_table"].obs["MF_niche"] = "None"
        gdf = sdata["transformed_cif_tiles"].loc[high_cif]


        for i in combined_polygons.index:
            poly = combined_polygons["geometry"][i]
            outcome = gdf.within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            sdata.tables["cif_512_table"].obs.loc[cell_ids, "MF_niche"] = i
            
        niches = combined_polygons
        cells = sdata["new_cell_boundaries"]

        cells_w_niches = geopandas.sjoin_nearest(cells, niches, distance_col="distances", how="left", max_distance=30000)
        
        cells_w_niches = cells_w_niches[["geometry", "index_right", "distances"]]


        #difference between two lists

        def diff(first, second):
            second = set(second)
            return [item for item in first if item not in second]

        #get duplicates in a list python

        def get_duplicates(lst):
            seen = set()
            seen_add = seen.add
            # adds all elements it doesn't know yet to seen and all other to seen_twice
            seen_twice = set( x for x in lst if x in seen or seen_add(x) )
            # turn the set into a list (as requested)
            return list( seen_twice )

        #get the indexes of the duplicates

        duplicates = get_duplicates(cells_w_niches.index)
        
        df_dropped = cells_w_niches.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        df_dropped = df_dropped.rename(columns={"index_right": "nearest_niche", "distances": "distance_to_niche"})
        df_dropped = df_dropped[["nearest_niche", "distance_to_niche"]]
        df = sdata.tables["table"].obs.merge(df_dropped, left_index=True, right_index=True)
        sdata.tables["table"].obs["nearest_niche"] = df["nearest_niche"]
        sdata.tables["table"].obs["distance_to_niche"] = df["distance_to_niche"]
        
        sdata.tables["table"].obs["nearest_niche"].fillna("None", inplace=True)
        sdata.tables["table"].obs["nearest_niche"] = sdata.tables["table"].obs["nearest_niche"].astype("str")
        sdata.tables["table"].obs["nearest_niche"] = sdata.tables["table"].obs["nearest_niche"].astype("category")
        
        niches = combined_polygons
        cells = sdata["transformed_cif_tiles"]

        cells_w_niches = geopandas.sjoin_nearest(cells, niches, distance_col="distances", how="left", max_distance=30000)
        cells_w_niches = cells_w_niches[["geometry", "index_right", "distances"]]

        def diff(first, second):
            second = set(second)
            return [item for item in first if item not in second]

        #get duplicates in a list python

        def get_duplicates(lst):
            seen = set()
            seen_add = seen.add
            # adds all elements it doesn't know yet to seen and all other to seen_twice
            seen_twice = set( x for x in lst if x in seen or seen_add(x) )
            # turn the set into a list (as requested)
            return list( seen_twice )

        #get the indexes of the duplicates

        duplicates = get_duplicates(cells_w_niches.index)

        df_dropped = cells_w_niches.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        df_dropped = df_dropped.rename(columns={"index_right": "nearest_niche", "distances": "distance_to_niche"})
        df_dropped = df_dropped[["nearest_niche", "distance_to_niche"]]

        df = sdata.tables["cif_512_table"].obs.merge(df_dropped, left_index=True, right_index=True)

        sdata.tables["cif_512_table"].obs["nearest_niche"] = df["nearest_niche"]
        sdata.tables["cif_512_table"].obs["distance_to_niche"] = df["distance_to_niche"]
        
        sdata.tables["cif_512_table"].obs["niche_status"] = "non-niche"
        #if distance to niche is 0 then it is a neighbour
        sdata.tables["cif_512_table"].obs.loc[sdata.tables["cif_512_table"].obs["distance_to_niche"] < 1, "niche_status"] = "neighbour"
        #if MF_niche is not none then it is a niche
        sdata.tables["cif_512_table"].obs.loc[sdata.tables["cif_512_table"].obs["MF_niche"] != "None", "niche_status"] = "niche"
        
        #make nearest niche a category and fill nan with none
        
        
        sdata.tables["cif_512_table"].obs["nearest_niche"].fillna("None", inplace=True)
        sdata.tables["cif_512_table"].obs["nearest_niche"] = sdata.tables["cif_512_table"].obs["nearest_niche"].astype("str")
        sdata.tables["cif_512_table"].obs["nearest_niche"] = sdata.tables["cif_512_table"].obs["nearest_niche"].astype("category")
        
        
        niches = sdata["transformed_bone"]
        cells = sdata["transformed_cif_tiles"]

        cells_w_niches = geopandas.sjoin_nearest(cells, niches, distance_col="distances", how="left", max_distance=30000)

        cells_w_niches = cells_w_niches[["geometry", "index_right", "distances"]]

        def diff(first, second):
            second = set(second)
            return [item for item in first if item not in second]

        #get duplicates in a list python

        def get_duplicates(lst):
            seen = set()
            seen_add = seen.add
            # adds all elements it doesn't know yet to seen and all other to seen_twice
            seen_twice = set( x for x in lst if x in seen or seen_add(x) )
            # turn the set into a list (as requested)
            return list( seen_twice )

        #get the indexes of the duplicates

        duplicates = get_duplicates(cells_w_niches.index)

        df_dropped = cells_w_niches.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")

        df_dropped = df_dropped.rename(columns={"index_right": "nearest_bone", "distances": "distance_to_bone"})

        df_dropped = df_dropped[["nearest_bone", "distance_to_bone"]]

        df = sdata.tables["cif_512_table"].obs.merge(df_dropped, left_index=True, right_index=True)

        sdata.tables["cif_512_table"].obs["nearest_bone"] = df["nearest_bone"]
        sdata.tables["cif_512_table"].obs["distance_to_bone"] = df["distance_to_bone"]
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_fibrotic_distance.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_fibrotic_distance.zarr"), overwrite=True)
        else:
            return sdata
    
    def tile_cellularity_512(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata

        for i in sdata["transformed_cif_tiles"].index:

            print(i)
            
            tile = sdata["transformed_cif_tiles"]["geometry"][i]

            #get the number of cells in the tile
            outcome = sdata["new_cell_boundaries"].within(tile)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            
            sdata.tables["cif_512_table"].obs.loc[i, "num_cells"] = len(cell_ids)
            
            #get the number of cells in the tile
            outcome = sdata["transformed_adipocytes_circles"].within(tile)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            
            sdata.tables["cif_512_table"].obs.loc[i, "num_adipocytes"] = len(cell_ids)
            
            #get the percentage area of fat in the tile
            
            result = geopandas.overlay(GeoDataFrame({"geometry": [tile]}), sdata["transformed_adipocytes"], how='intersection')
            
            per_fat = (result.area.sum() / tile.area ) * 100
            
            sdata.tables["cif_512_table"].obs.loc[i, "per_fat"] = per_fat
            
            #get the percentage of cells in the tile
            
            #theoretically I should use new_cell_boundaries with_fill_adipocytes and remove the adipocytes 
            
            result = geopandas.overlay(GeoDataFrame({"geometry": [tile]}), sdata["new_cell_boundaries"], how='intersection')
            
            per_cells = (result.area.sum() / tile.area ) * 100
            
            sdata.tables["cif_512_table"].obs.loc[i, "per_cells"] = per_cells
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_512_cellularity.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_512_cellularity.zarr"), overwrite=True)
        else:
            return sdata
    
    def cellularity_it(self, zarr_path=None, sdata=None, save=False):
        
        if zarr_path != None:
            sdata = self.initialise_sdata_from_zarr(zarr_path)
        else:
            sdata = sdata
            
        for i in sdata["intertrabecular_regions"].index:

            print(i)
            
            tile = sdata["intertrabecular_regions"]["geometry"][i]

            #get the number of cells in the tile
            outcome = sdata["new_cell_boundaries"].within(tile)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            
            sdata.tables["intertrabecular_regions_table"].obs.loc[i, "num_cells"] = len(cell_ids)
            
            #get the number of cells in the tile
            outcome = sdata["transformed_adipocytes_circles"].within(tile)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            
            sdata.tables["intertrabecular_regions_table"].obs.loc[i, "num_adipocytes"] = len(cell_ids)
            
            #get the percentage area of fat in the tile
            
            result = geopandas.overlay(GeoDataFrame({"geometry": [tile]}), sdata["transformed_adipocytes"], how='intersection')
            
            per_fat = (result.area.sum() / tile.area ) * 100
            
            sdata.tables["intertrabecular_regions_table"].obs.loc[i, "per_fat"] = per_fat
            
            #get the percentage of cells in the tile
            
            #theoretically I should use new_cell_boundaries with_fill_adipocytes and remove the adipocytes 
            
            result = geopandas.overlay(GeoDataFrame({"geometry": [tile]}), sdata["new_cell_boundaries"], how='intersection')
            
            per_cells = (result.area.sum() / tile.area ) * 100
            
            sdata.tables["intertrabecular_regions_table"].obs.loc[i, "per_cells"] = per_cells
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_it_cellularity.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_it_cellularity.zarr"), overwrite=True)
        else:
            return sdata    
        
        
        
    
    
        
        
            
        
        
        
        
    
            
            

                
            
            
            
            
        
            
        
        
        
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
 
    #add in negative polygons
    
    
    #add in intratrebecular polygons
    
    
    #mark cells to remove under negative polygons
    

    #mark cells in each intratrebecaular region
    
    
    #add in cif scores
    
    
    
    
    
    
    
    
        
    # def add_meg_annotations(self, sdata):
    #     meg = pd.read_csv(self.meg_path)
    #     meg_pheno = pd.read_csv(self.meg_pheno_path)
        
    #     sdata.add_table("meg", TableModel(meg))
    #     sdata.add_table("meg_pheno", TableModel(meg_pheno))
        
    #     return sdata
    
    # def add_he_image(self, sdata):
    #     sdata.add_image("he_image", spatialdata_io.image.imread(self.he_image_path))
        
    #     return sdata
    
    # def add_annotations(self, sdata):
    #     annotation = pd.read_csv(self.annotation_path + self.sample_key + ".csv")
    #     sdata.add_table("annotations", TableModel(annotation))
        
    #     return sdata
    
    # def add_cif(self, sdata):
    #     cif = pd.read_csv(self.cif_path + self.sample_key + ".csv")
    #     sdata.add_table("cif", TableModel(cif))
        
    #     return sdata
    
    # def add_neg_annotations(self, sdata):
    #     neg_annotation = pd.read_csv(self.neg_annotation_save + self.sample_key + ".csv")
    #     sdata.add_table("neg_annotations", TableModel(neg_annotation))
        
    #     return sdata
    
    # def add_tissue_poly(self, sdata):
    #     tissue_poly = pd.read_csv(self.tissue_poly)
    #     sdata.add_table("tissue_poly", TableModel(tissue_poly))
        
    #     return sdata
    
    # def add_all_data(self, sdata):
    #     sdata = self.add_meg_annotations(sdata)
    #     sdata =
            
        

