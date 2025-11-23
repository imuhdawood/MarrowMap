# import dask
# from pathlib import Path
# dask.config.set({'dataframe.query-planning': False})

# import sopa
# from scipy.sparse import csr_matrix
# #import spatialdata_xenium_explorer
# from sopa.utils.data import uniform
# import sopa.segmentation
# from spatialdata_io import xenium
# import spatialdata as sd
# import spatialdata_plot
# import os
# import matplotlib.pyplot as plt
# import json
# from shapely import MultiPolygon, Polygon, Point, make_valid
# from spatialdata.models import ShapesModel, TableModel
# import geopandas as gpd
# from spatialdata.transformations import (
#     Identity,
#     get_transformation_between_coordinate_systems,
#     set_transformation,
# )
# from spatialdata import transform
# import pandas as pd

# import numpy as np
# import json

# from geopandas import GeoDataFrame
# from shapely.geometry import Polygon
# import shapely

# from shapely import LineString, Point
# import shapely.wkt
# from shapely import get_coordinates

# import anndata
# import scanpy as sc

# from spatialdata.models import TableModel
# import matplotlib.pyplot as plt

# import numpy as np
# from spatialdata.models import ShapesModel

# import math
# import json
# import geopandas
# import anndata as ad

# from shapely import geometry 



# Standard library
from pathlib import Path
import json
import math
import os

# Array/data handling
import dask
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# Spatial / geometry tooling
import geopandas as gpd
from geopandas import GeoDataFrame
import shapely
import shapely.wkt
from shapely import LineString, MultiPolygon, Point, Polygon, geometry, get_coordinates, make_valid
from shapely.ops import unary_union
import sopa
import sopa.segmentation
from sopa.utils.data import uniform
import spatialdata as sd
from spatialdata import transform
from spatialdata.models import ShapesModel, TableModel
from spatialdata.transformations import (
    Identity,
    get_transformation_between_coordinate_systems,
    set_transformation
)

# Single-cell / AnnData helpers
import anndata as ad
import scanpy as sc


import sys
sys.path.append('/well/rittscher/users/qwi813/xenium_paper')
from configuration import STConfig
from wrappers.sdata_initalizer import SDataInitalizer
from utils.utils_sdata import get_intrinsic_cs, to_intrinsic
from utils.utils_geom import clean_geometry, rescale_polygons


class SDataCustomizer(SDataInitalizer):
    def __init__(self, config: STConfig,
                 row_idx=0,
                 sample_key='10693_R2',
                 key_col='sample_key',
                 verbose = True
                 ):
        
        super().__init__(config, row_idx=row_idx, sample_key=sample_key, key_col=key_col)
        self.pth_landmark_matrix = config.pth_landmark_matrix
        self.verbose = verbose
        self.MPP = config.MPP
        self.ADIPO_BUFFER_RADIUS = config.ADIPO_BUFFER_RADIUS
        self.CELL_BONE_MAX_DISTANCE = config.CELL_BONE_MAX_DISTANCE
        self.PERITRABECULAR_DISTANCE_THRESHOLD = config.PERITRABECULAR_DISTANCE_THRESHOLD
        self.ENDOSTEAL_DISTANCE_THRESHOLD = config.ENDOSTEAL_DISTANCE_THRESHOLD

        self.pth_pos_annotation_save = config.pth_pos_annotation_save
        self.pth_neg_annotation_save = config.pth_neg_annotation_save
        self.pth_cell_annotations = config.pth_cell_annotations # file containing cell annotations
        self.pth_feature_map = config.pth_feature_map
        self.pth_bone_segmentation = config.pth_bone_segmentation
        self.pth_adipocytes_segmentation = config.pth_adipocytes_segmentation
    
    def _logs(self, message: str) -> None:
        if self.verbose:
            print(message)

    def _resolve_sdata_inputs(self, zarr_path=None, sdata=None):
        if zarr_path:
            return self.initialise_sdata_from_zarr(zarr_path), zarr_path
        if sdata is not None:
            return sdata, None
        default_path = os.path.join(self.save_folder, f"{self.sample_key}.zarr")
        self._logs(f'Both zarr_path and sdata is None \n loading sample {default_path}')
        return self.initialise_sdata_from_zarr(default_path), default_path
    

    def add_he_image(
            self, zarr_path: str | None =None
            ):

        """
        Generate the sopa CLI command to attach an aligned H&E image to an sdata zarr
        using the instance's sample_key and landmark-matrix naming convention.
        """
        # adjust name
        l_he_key = self.he_key.replace("x", "X")
        if not zarr_path:
            zarr_path = os.path.join(self.save_folder, self.sample_key + '.zarr')
        
        bash_command = f"sopa explorer add-aligned {zarr_path} {self.he_image_path} {self.pth_landmark_matrix}/{l_he_key}_{self.sample_key}_matrix.csv --original-image-key morphology_focus"
        if self.verbose:
            print(bash_command)
        return bash_command
    
    
    def add_megakaryocytes(
            self,
            zarr_path = None,
            sdata = None,
            save=False,
            pattern = r'fromx_(\d+)_fromy_(\d+)_tox_(\d+)_toy_(\d+)\.png$',
            postfix = 'with_megs'
            ):
        
        """
        Add megakaryocyte shapes to the sample’s SpatialData.

        Loads the zarr, reads meg polygons (GeoJSON) and phenotype CSV, extracts
        bbox centroids from filenames (via `pattern`), assigns phenotypes, joins
        points to polygons, transforms to 'global' from H&E (self.he_key), stores
        as `sdata.shapes["transformed_megs"]`, and optionally saves.

        Args:
            save (bool): Write updated zarr if True.
            pattern (str): Regex to parse xmin,ymin,xmax,ymax from filenames.
            postfix (str): Suffix to append before ".zarr" when saving.

        Returns:
            spatialdata.SpatialData: Updated object with 'transformed_megs'.
        """

        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)   
        # Load polygons of each meg
        with open(self.meg_path) as f:
            data = json.load(f)
        megs_poly = [Polygon(f['geometry']['coordinates'][0]) for f in data['features']]

        # extracting meg cells centroids
        meg_pheno = pd.read_csv(self.meg_pheno_path)
        coords = np.array(meg_pheno['name'].str.extract(pattern).astype(int))
        Cx = coords[:,[0, 2]].mean(1) # mmin and mmax
        Cy = coords[:,[1, 3]].mean(1) # ymin and ymax
        meg_centroids = list(zip(Cx.tolist(), Cy.tolist()))
        meg_pheno['centroid'] = meg_centroids
        meg_points = [Point(xy) for xy in meg_centroids]
        points_gdf = gpd.GeoDataFrame({"geometry": meg_points, "phenotype": meg_pheno['membership']})
        
        # spatial join
        gdf = gpd.GeoDataFrame({"geometry": megs_poly})
        points_within = gpd.sjoin(gdf, points_gdf, how="left", predicate='contains')
        points_within.drop(columns=["index_right"], inplace=True)
        points_within.drop_duplicates(subset="geometry", inplace=True)

        # apply coordinates transformation
        new_gdf = ShapesModel.parse(points_within, transformations={"global": Identity()})
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        transformed_gdf = transform(new_gdf, affine_he_to_xen, "global")
        sdata.shapes["transformed_megs"] = transformed_gdf
        set_transformation(sdata.shapes["transformed_megs"], Identity(), to_coordinate_system="global")
        
        if save == True:
            if zarr_path:
                zarr_path = os.path.join(self.save_folder, self.sample_key + '.zarr')
            zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")    
            sdata.write(zarr_path, overwrite=True)
            self._logs(f"File save {zarr_path}")
        return sdata
    
    
    def add_megs_to_segmentation(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        postfix = 'with_new_seg'
        ):

        """
        Create an updated cell segmentation that integrates HE-derived megakaryocytes.
        Validates and rescales cell/nucleus boundaries to pixel space, transforms megakaryocytes
        into the same frame, subtracts megakaryocytes from cells, resolves resulting multipolygons
        using nuclei (keeping full cells, merging small fragments to nearest nuclei), flags artefacts,
        and assigns a `cell_status` per polygon. Writes the result to
        `sdata.shapes['new_cell_boundaries']`, saves a CSV of statuses, and optionally writes a
        `*_with_new_seg.zarr`.

        Args:
            zarr_path: Optional path to the source Zarr; used when saving.
            sdata: Optional SpatialData object; resolved if not provided.
            save: If True, writes out the updated Zarr; otherwise returns updated `sdata`.

        Returns:
            The updated `sdata` when `save=False`.
        """
        
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)

        self._logs("Validating megakaryocyte polygons")
        if "transformed_megs" not in sdata.shapes:
            raise KeyError("Add megakaryocytes first by calling add_megakaryocytes() methods")
        
        meg_gdf = clean_geometry(
            sdata.shapes["transformed_megs"],
            self._logs
            )
        sdata.shapes["transformed_megs"]  = meg_gdf
        
        self._logs("Validating and transforming cell boundaries")
        cell_gdf = clean_geometry(
            sdata.shapes["cell_boundaries"],
            self._logs
            )
        
        # Convert xenium MPP to pixel space
        cell_gdf = rescale_polygons(cell_gdf, self.MPP)
        set_transformation(cell_gdf, Identity(), "global")
        sdata.shapes["cell_boundaries"] = cell_gdf
        
        self._logs("Validating and transforming cell boundaries")

        #geo_df into same space as cells
        megs_boundaries = to_intrinsic(sdata, sdata["transformed_megs"], cell_gdf)
        cell_boundaries = cell_gdf.copy()  
        cell_boundaries["cell_id"] = cell_boundaries.index

        gdf_join = cell_boundaries.overlay(megs_boundaries, how="difference")       
        #check for multipolygons in gdf_join - these will be left over cells from the overlaps
        mp = gdf_join.loc[gdf_join.geometry.geom_type != "Polygon"]
        
        self._logs(f"Left over polygon after HE megs and Xenium cell boundaries overlap {len(mp)}")
        #explode each multipolygon into it's list of polygons
        mp = mp.explode().reset_index(drop=False)
        #make a new column by append the index to the cell_id
        mp["new_cell_id"] = mp["cell_id"] + "_" + mp["level_1"].astype(str)
        #make new_cell_id the index
        mp.index = list(mp["new_cell_id"])
        #drop some columns from mp
        mp = mp.drop(columns=["level_0", "level_1", "new_cell_id"])
        invalid = mp[~mp.geometry.is_valid]
        self._logs(f'Cell with invalid polygon boundary # {len(invalid)}')
        
        self._logs("Check if polygons contain nuclei")
        
        nuclei = rescale_polygons(
            sdata.shapes["nucleus_boundaries"], self.MPP
        )
        
        #now change teh transformation to identity because we have changed them manually
        set_transformation(nuclei, Identity(), "global")
        
        #add these back into sdata
        sdata["nucleus_boundaries"] = nuclei
        
        #get the centroid of each polygon
        nuc_gdf = gpd.GeoDataFrame(
            {'radius':1},
            geometry=nuclei.geometry.centroid,
            index=nuclei.index
            )
        
        #remove multipolygons from new cell segmentation
        poly_gdf_join = gdf_join.loc[gdf_join.geometry.geom_type == "Polygon"]
        #add full cells to poly_gdf_join these are not altered
        # Both original and HE based model boundaries overlaps
        poly_gdf_join["cell_status"] = "original"
        
        #get nuclei which are in megakaryocytes
        meg_gdf["nuclei"] = None

        for i in meg_gdf["geometry"].index:
            pol = meg_gdf["geometry"][i]
            #pol_gdf = gpd.GeoDataFrame({"geometry": [pol]})
            nuc_meg = nuc_gdf.within(pol)
            #get indexes where nuc_meg is true
            idx = list(nuc_meg[nuc_meg == True].index)
            meg_gdf["nuclei"][i] = idx
        #now remove the nuclei which are in megakaryocytes from nuc_gdf
        nuc_in_megs = meg_gdf['nuclei'].explode().dropna().unique()
        sub_nuc_gdf = nuc_gdf.drop(nuc_in_megs)
        
        # if the MultiPolygon exploded into cells contains a nucleus,
        # put it back into cell_boundaries with a note to say altered
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
        
        #check that all the full cells are unique
        self._logs(f"Are full cells unique:, {full_cells.shape[0] == full_cells['cell_id'].nunique()}")
        
        #drop full_cell column from full_cells
        full_cells = full_cells.drop(columns=["full_cell"])
        #add column
        full_cells["cell_status"] = "altered"
        #concatenate full_cells and poly_gdf_join together
        full_cells_and_poly_gdf_join = pd.concat([poly_gdf_join, full_cells])

        self._logs("Deal with small cells and artefacts")
        # Step 1: selecting cells that are not full cells
        artefacts_gdf = mp[mp["full_cell"] == False]
        # Step 2: Cells with area greater that median area of all cells are thrown away
        median_cell_area = np.median(cell_gdf.area)/2
        art_area = artefacts_gdf.area
        small_cells = artefacts_gdf[art_area < median_cell_area]
        artefacts = artefacts_gdf[art_area >= median_cell_area]
        artefacts.drop(columns=["full_cell"], inplace=True)
        artefacts["cell_id"] = artefacts.index.tolist()
        artefacts["cell_status"] = "artefact"

        # Step 3: Find the nearest nucleus, that is not in a megakaryocyte, to each small cell
        # merge small cell into the cell that contains that nucleus

        #not sure which point in polygon uses to calculate distance
        #relying on nuclei having same names as cells here
        nn_join = small_cells.sjoin_nearest(sub_nuc_gdf, max_distance=100, distance_col="distances")
        
        #okay now merge small cells into the cells in index_right
        cells_for_merge = list(set(nn_join["index_right"]))
        
        #check if each nuclei in cells_for_merge has a equivalent cell in index of full_cells_and_poly_gdf_join
        #check if list contains all elemnets of another list
        set(cells_for_merge).issubset(list(full_cells_and_poly_gdf_join["cell_id"]))
        
        self._logs("Conglomerating all polygons together into new segmentation")
        
        small_cells_full_cells_and_poly_gdf_join = full_cells_and_poly_gdf_join.copy()

        for i in cells_for_merge:
            #get the polygons where index_right = i
            #get the small cells out 
            nn_sub = nn_join[nn_join["index_right"] == i]
            nn_sub = nn_sub.drop(columns=["distances", "radius", "index_right", "full_cell"])
            #print(nn_sub)
            
            #get the large cell taht it needs to be merged into
            nn_large = full_cells_and_poly_gdf_join[full_cells_and_poly_gdf_join["cell_id"] == i]
            #print(nn_large)
            
            #concatenate these two polygons into one geodataframe
            merged = pd.concat([nn_sub, nn_large])

            #now merge these polygons into one polygon using unary_union
            #expand so we can merge then de-expand - otherwise end up with a multipolygon
            union_poly = merged.buffer(10).unary_union.buffer(-10)
            #print(union_poly)
            
            small_cells_full_cells_and_poly_gdf_join.index = list(small_cells_full_cells_and_poly_gdf_join["cell_id"])
            small_cells_full_cells_and_poly_gdf_join.loc[i, "geometry"] = union_poly
            small_cells_full_cells_and_poly_gdf_join.loc[i,"cell_status"]= "altered_with_small_cells"
            
        self._logs("add in artefacts and megakaryocytes")
        
        #okay now I just need to add in the artefact cells and megakaryocytes
        meg_gdf["cell_status"] = "megakaryocyte"
        meg_gdf.drop(columns="nuclei", inplace=True)
        meg_gdf["cell_id"] = list(meg_gdf.index)
        
        #now concatenate all together
        new_cell_gdf = pd.concat([small_cells_full_cells_and_poly_gdf_join, artefacts, meg_gdf])
        
        #if cell_status isn't original append new to cell_id and make index - else leave it so doesn't have to be updated
        for i in list(new_cell_gdf.index):
            if new_cell_gdf["cell_status"][i] != "original":
                new_cell_gdf.loc[i, "cell_id"] = str(i) + "_new"
                
        new_cell_gdf.index = list(new_cell_gdf["cell_id"])
        
        invalid = new_cell_gdf.loc[~new_cell_gdf.geometry.is_valid]
        
        #parse this with identity matrix
        set_transformation(new_cell_gdf, Identity(), "global")
        sdata.shapes["new_cell_boundaries"] = new_cell_gdf
        sdata["new_cell_boundaries"].assign(cell_status = list(new_cell_gdf["cell_status"]))
        
        cell_meta_file = Path(self.save_folder) / f"cells_meta_{self.sample_key}.csv"
        new_cell_gdf.to_csv(cell_meta_file, index=False)
        self._logs(f'Updated cell geometry with meta data saved \n {cell_meta_file}')

        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
    
    
    def aggregate_with_sopa(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        postfix = 'sopa_aggregated'
        ):

        """
        Run Sopa aggregation on `new_cell_boundaries`, clean multipolygons, then update
        the transcript table and cell metadata before optionally writing the augmented Zarr.
        """
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)
        print(sdata)  
        #sdata.path = None # setting the out path to none so that the file is not overwritten
        if "new_cell_boundaries" not in sdata.shapes:
            raise KeyError("new_cell_boundaries not present in sdata")
        
        #extract the largest polygon from each multipolygon
        cell_gdf = sdata["new_cell_boundaries"]
        multi_poly_cells = cell_gdf[cell_gdf.geometry.geom_type=='MultiPolygon'].index.tolist()
        if len(multi_poly_cells) > 0:
            for cell_id in multi_poly_cells:
                if cell_gdf["geometry"][cell_id].geom_type == "MultiPolygon":
                    #choose the biggest polygon from each 
                    largest_poly = max(list(cell_gdf["geometry"][cell_id].geoms), key=lambda a: a.area)
                    cell_gdf.geometry[cell_id] = largest_poly
            
            #  Updating geometry in sdata
            sdata["new_cell_boundaries"] = cell_gdf

        # Loading the cell meta to update andata as well and add a Flag
        cell_meta_file = Path(self.save_folder) / f"cells_meta_{self.sample_key}.csv"
        new_cell_gdf = pd.read_csv(cell_meta_file)
        self._logs(f"# of Cells {new_cell_gdf.shape}")
        
        old_anndata = sdata.tables["table"].copy()
        # sopa throw exception and need worked after dropping cell_id 
        #https://gustaveroussy.github.io/sopa/tutorials/api_usage/
        sdata.path = None # setting the out path to none so that the file is not overwritten
        aggregator = sopa.segmentation.Aggregator(sdata, image_key="morphology_focus", shapes_key="new_cell_boundaries", overwrite=False)
        #need to stop it changing the indices
        aggregator.compute_table(gene_column="feature_name", average_intensities=False)

        #change index (know this is fine because area the same)
        sdata["new_cell_boundaries"].index = new_cell_gdf["cell_id"].tolist()

        table = sdata.tables["table"].copy()
        table.obs.index = new_cell_gdf["cell_id"].tolist()
        table.obs.loc[:, ["cell_status", "cell_id"]] = new_cell_gdf[["cell_status", "cell_id"]].to_numpy()
        table = table[new_cell_gdf["cell_id"], old_anndata.var_names]
        table.obs["region"] = pd.Categorical(["new_cell_boundaries"] * table.n_obs,
                                            categories=["new_cell_boundaries"])
        sdata.tables["table"] = table
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
        
    #add cif scores
    
    def add_cif_scores(
        self,
        zarr_path=None,
        window_size=512,
        sdata=None,
        save=False,
        cif_path=None,
        save_transformed_cif = False,
        postfix = 'with_512_cif'
        ):

        """
        Load CIF tile scores, keep those whose footprint lies within `tissue_poly`,
        convert them into a SpatialData `Shapes` element, transform to the Xenium
        global coordinates, optionally persist both the transformed CSV and an updated Zarr.
        """

        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)
        
        if cif_path != None:
            self.cif_path = cif_path
        
        cif_score = pd.read_csv(os.path.join(self.cif_path, self.he_key + ".txt"), delim_whitespace = True)

        # Selecting tiles within tissue polygon of this sample as multiple samples are 
        # scanned on a single H&E slides and tissue polygon contains polygon for this sample
        P = shapely.wkt.loads(self.tissue_poly)
        points = gpd.GeoSeries(gpd.points_from_xy(cif_score["X"], cif_score["Y"]))
        mask = points.within(P)
        sub_cif_score = cif_score[mask]
         
        tile_list = []  

        for i in list(sub_cif_score.index):
            minx = sub_cif_score["X"][i]
            miny = sub_cif_score["Y"][i]
            maxx = minx + window_size
            maxy = miny + window_size
            bbox_tuple = [minx, miny, maxx, maxy]
            bbox_polygon = geometry.box(*bbox_tuple)
            tile_list.append(bbox_polygon)
        
        cif_gdf = gpd.GeoDataFrame({"geometry": tile_list, "cif_score": sub_cif_score['Score']})
        cif_gdf = ShapesModel.parse(cif_gdf, transformations={"global": Identity()})
        
        # transformation matrix between sdata and H&E image to apply to cif_gdf
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        transformed_cif_gdf = transform(cif_gdf, affine_he_to_xen, "global")
        sdata.shapes["transformed_cif_tiles"] = transformed_cif_gdf
        # updating metadata
        set_transformation(sdata.shapes["transformed_cif_tiles"], Identity(), to_coordinate_system="global")
        
        if save_transformed_cif != False:
            transformed_cif_gdf.to_csv(os.path.join(save_transformed_cif, self.sample_key + "_512_transformed.csv"), index=False)
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
    
    
    def add_meg_phenotypes(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        postfix = 'with_meg_pheno'
        
        ):

        """
        Assuming that megakaryocyte shapes are alreadding added to SpatialData object using
        `add_megakaryocyte` function. The function copy meg phenotype labels
        into `sdata.tables["table"].obs['meg_phenotype']`, marking all
        other cells as `non_meg`, and optionally write the updated Zarr store.
        """
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)
        
        # Adding meg phenotypes to adata.obs table by matching indices of meg_df
        # with sdata.tables["table"].obs where they are tagged as meg-idx_new
        meg_gdf = sdata["transformed_megs"]
        meg_cell_ids = [f'{i}_new' for i in meg_gdf.index.tolist()]
        meg_gdf["cell_id"] = meg_cell_ids
        meg_gdf.index = meg_cell_ids
        
        sdata.tables["table"].obs["meg_phenotype"] = "non_meg"
        sdata.tables['table'].obs.loc[meg_cell_ids,'meg_phenotype'] = meg_gdf['phenotype'].tolist()
            
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
    
    def load_manual_annotations(
        self,
        csv_subscript="_N_C.csv"
        ):

        """
        Read Xenium Explorer exported manual-annotation CSVs (coordinates already in microns),
        convert them back to pixel space using`self.MPP` as sampling factor, and return a
        GeoDataFrame containing each selection polygon along with its region ID.
        """

        ann_path = os.path.join(self.annotation_path, self.sample_key)
        selections = pd.read_csv(os.path.join(ann_path, self.sample_key + csv_subscript), sep=",", comment='#')

        # Make index a column
        selections = selections.reset_index()
        # Make the first row the header
        selections.columns = selections.iloc[0]
        selections = selections[1:]
        selections[["X", "Y"]] = selections[["X", "Y"]].astype(float).div(self.MPP)

        n_selections = selections["Selection"].nunique()
        regions = []
        region_ids = []
        for i in range(n_selections):
            sel = selections[selections["Selection"] == selections["Selection"].unique()[i]]
            poly = Polygon(sel[["X", "Y"]].values.tolist())
            regions.append(poly)
            # Get regions id and add it to a list
            num = selections["Selection"].unique()[i].split(" ")[-1]
            region_ids.append(num)
        gdf = gpd.GeoDataFrame({"geometry": regions, "region_id": region_ids})
        gdf = ShapesModel.parse(gdf, transformations={"global": Identity()})
        return gdf

    def add_negative_annotations(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        save_neg_annotations=None,
        postfix = 'with_neg_poly'
        ):

        """
        Load the manual “negative” annotations, store them as the `crushed_regions`
        shapes layer, export their H&E-space coordinates if requested, flag any cells
        whose boundaries fall inside those polygons as `annotation = "remove"`, then
        optionally persist the updated SpatialData Zarr.
        """

        
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)

        if save_neg_annotations == None:
            save_neg_annotations = self.pth_neg_annotation_save
        
        gdf = self.load_manual_annotations(csv_subscript="_N_C.csv")

        sdata.shapes["crushed_regions"] = gdf

        affine_glob_to_he = get_transformation_between_coordinate_systems(sdata, "global", sdata.images[self.he_key])
        transformed_gdf = transform(gdf, affine_glob_to_he, "he_space")
        
        if save_neg_annotations:
            os.makedirs(save_neg_annotations, exist_ok=True)
            transformed_gdf.to_csv(os.path.join(save_neg_annotations, self.sample_key + "_N.csv"), index=False)
        
        # now add in negative annotations to the cell table
        sdata.tables["table"].obs["annotation"] = "keep"

        cells_to_remove = []

        for i in list(sdata["crushed_regions"].index):
            poly = sdata["crushed_regions"].geometry[i]
            outcome = sdata["new_cell_boundaries"].within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            cells_to_remove.extend(cell_ids)
        
        sdata.tables["table"].obs.loc[cells_to_remove,"annotation"] = "remove"
        self._logs(f"Cells to remove based on negative annotations: {len(cells_to_remove)}")

        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
        
    def add_positive_annotations(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        save_pos_annotations=None,
        postfix = 'with_pos_poly'
        ):
        
        """
        Load the “positive” manual annotations (intertrabecular regions), add them as a
        shapes layer, optionally export their H&E-space coordinates, tag cells whose
        boundaries fall inside those polygons with the corresponding region ID (others
        stay `non_intertrabecular`), log the counts, and optionally write the Zarr.
        """

        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)

        if save_pos_annotations == None:
            save_pos_annotations = self.pth_pos_annotation_save 

        gdf = self.load_manual_annotations(csv_subscript="_P_C.csv")  
        sdata.shapes["intertrabecular_regions"] = gdf
        
        #get the he transformed version to save for Hosuk

        affine_glob_to_he = get_transformation_between_coordinate_systems(sdata, "global", sdata.images[self.he_key])
        transformed_gdf = transform(gdf, affine_glob_to_he, "he_space")
        
        if save:
            os.makedirs(save_pos_annotations, exist_ok=True)
            transformed_gdf.to_csv(os.path.join(save_pos_annotations, self.sample_key + "_P.csv"), index=False)

        # Adding in positive intratrebecular annotations

        sdata.tables["table"].obs["it_regions"] = "non_intertrabecular"

        it_cells = []
        cell_ids_it = []

        for i in list(sdata["intertrabecular_regions"].index):
            poly = sdata["intertrabecular_regions"].geometry[i]
            outcome = sdata["new_cell_boundaries"].within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            cell_ids_it.extend(cell_ids)
            it = [i]*len(cell_ids)
            it_cells.extend(it) 
        
        df = pd.DataFrame(it_cells, columns=["it_regions"], index = cell_ids_it)
        df.index = cell_ids_it
        
        df2 = df[~df.index.duplicated(keep='first')]
        sdata.tables["table"].obs.update(df2)
        
        self._logs(f'# of ITS regions {sdata.tables["table"].obs["it_regions"].nunique()}')

        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata 
        
    
    def add_metadata(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        postfix = 'with_metadata'
        ):

        """
        Normalize cell IDs (strip `_new`, append sample key), propagate those IDs to
        the `new_cell_boundaries` index, then attach run/sample metadata
        (`mutation_status`, `run`, `study_id`, diagnosis fields, meg phenotypes/ ITS labels) to
        every row in `sdata.tables["table"].obs`. This duplicates cohort metadata per
        cell—inefficient but convenient when cell-level access is required—and can
        optionally persist the enriched Zarr.
        """

        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)

        sdata.tables["table"].obs["cell_id"] = sdata.tables["table"].obs["cell_id"].str.replace("_new", "")
        sdata.tables["table"].obs["cell_id"] = sdata.tables["table"].obs["cell_id"] + "_" + self.sample_key
        sdata["new_cell_boundaries"].index = list(sdata.tables["table"].obs["cell_id"])
        sdata.tables["table"].obs.index = list(sdata.tables["table"].obs["cell_id"])

        # Add other information. Through this is not a memory efficient way as for each cell these would be duplicated
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
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata

    
    def add_cell_annotations(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        anno_path = None,
        original_cluster_name="SCT_snn_res.0.3"
        ):

        """
        Import legacy cell-type annotations done at the start of the poject for EHA Conference,
        Adding those legacy annotations to the sdata object

        Note: These annotations are not used in current paper analyses.
        """

        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)

        if anno_path == None:
            anno_path = self.pth_cell_annotations
        
        if not os.path.exists(anno_path):
            self._logs(f"Annotation path does not exist: {anno_path}")
            return
        else:
            self._logs(f"Loading annotations from: {anno_path}")
            anno = pd.read_csv(anno_path)

        # make sample id column by extracting last 8 characters from cell_id
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
        
    def change_cell_annotations(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        map_path = None,
        original_cluster_name="SCT_snn_res.0.3",
        postfix = 'with_changed_annotations'
        ):

        """
        Import legacy cell-type annotations done at  round 2 of the poject. The approach used again was 
        only unsupervised cell type assignment. Adding that to sdata obejct for reference.
        
        Note: These annotations are not used in current paper analyses.
        """
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)

        if map_path == None:
            map_path = self.pth_feature_map
        
        if not os.path.exists(map_path):
            self._logs(f"Annotation path does not exist: {map_path}")
            return
        else:
            self._logs(f"Loading annotations from: {map_path}")       
            df_map = pd.read_csv(map_path)
            
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
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata 
    
    def add_bone_segmentation(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        postfix = 'with_bone'
        ):

        """
        Load the per-sample bone segmentation GeoJSON, convert it into a SpatialData
        Shapes element, transform it from the H&E coordinates into the Xenium global coords,
        add it under `sdata.shapes["transformed_bone"]`, and optionally persist
        an updated Zarr (suffix `_with_bone` by default).
        """

        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)
        
        pth_bone_segmentation = self.pth_bone_segmentation
        bone_geo_path = os.path.join(pth_bone_segmentation, self.sample_key + ".geojson")

        if not os.path.isfile(bone_geo_path):
            ValueError(f"No bone segmentation found for sample {self.sample_key} at {bone_geo_path}")

        with open(bone_geo_path) as f:
            dd = json.load(f)
            
        gdf = GeoDataFrame.from_features(dd)
        gdf = ShapesModel.parse(gdf, transformations={"global": Identity()})
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        transformed_gdf = transform(gdf, affine_he_to_xen, "global")
        
        transformed_gdf.index = [f"{self.sample_key}_B{i}" for i in range(len(transformed_gdf))]
        sdata.shapes["transformed_bone"] = transformed_gdf
        set_transformation(sdata.shapes["transformed_bone"], Identity(), to_coordinate_system="global")
            
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
    
    
    def add_adipocytes(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        postfix = 'with_unaggregated_adipocytes'
        ):

        """
        Load per-sample adipocyte polygons from GeoJSON, drop spurious LineStrings,
        register them as a Shapes layer, transform from the H&E frame into the Xenium
        global frame, and produce both full polygons and buffered circle surrogates.
        The circle layer is used to subtract adipocytes from `new_cell_boundaries`
        via a symmetric difference, then the remaining cells are recombined with either
        the circular surrogates or full polygons to produce three updated boundaries
        layers (`…without_adipocytes`, `…with_adipocytes`, `…with_full_adipocytes`).
        Optionally writes the augmented Zarr to disk.
        """

        
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)
        
        pth_adipocytes_segmentation = self.pth_adipocytes_segmentation
        adipo_geo_path = os.path.join(pth_adipocytes_segmentation, self.sample_key + ".geojson")

        if not os.path.isfile(adipo_geo_path):
            ValueError(f"No bone segmentation found for sample {self.sample_key} at {adipo_geo_path}")
        
        with open(adipo_geo_path) as f:
            data = json.load(f)
        
        geodf = GeoDataFrame.from_features(data["features"], columns=["geometry"])
        
        # Excluding this as these are less likely to be adipocytes
        geodf = geodf[geodf["geometry"].geom_type != "LineString"]
        # adding adipcyte geometry to sdata object as a shape layer 
        geodf.reset_index(drop=True, inplace=True)
        geodf = ShapesModel.parse(geodf, transformations={"global": Identity()})
        sdata.shapes["adipocytes"] = geodf
        
        affine_he_to_xen = get_transformation_between_coordinate_systems(sdata, sdata.images[self.he_key], "global")
        
        transformed_gdf = transform(geodf, affine_he_to_xen, "global")
        new_transformed_representation = transformed_gdf.copy()
        transformed_gdf.index = [f"{self.sample_key}_A{i}" for i in range(len(transformed_gdf))]
        sdata.shapes["transformed_adipocytes"] = transformed_gdf

        #new_transformed_gdf["geometry"] = transformed_gdf["geometry"].apply(lambda x: x.buffer(-50) if x.area > 8000 else x.buffer(-30))
        
        #now make each of those centers into a circle of radius 15 pixels and diameter 30 pixels
        centers = new_transformed_representation["geometry"].centroid
        new_transformed_representation["geometry"] = [Point(i) for i in centers]
        new_transformed_representation["geometry"] = new_transformed_representation["geometry"].apply(lambda x: x.buffer(self.ADIPO_BUFFER_RADIUS))
        new_transformed_representation.index = [f"{self.sample_key}_A{i}" for i in range(len(new_transformed_representation))]
        sdata.shapes["transformed_adipocytes_circles"] = new_transformed_representation
        # set transformations to identity as already in global coords
        set_transformation(sdata.shapes["transformed_adipocytes"], Identity(), to_coordinate_system="global")
        set_transformation(sdata.shapes["transformed_adipocytes_circles"], Identity(), to_coordinate_system="global")
        
        # Do not know the exact motivation but one reason could be to improve cell segmentation need to check with @Emily
        df1 = sdata["new_cell_boundaries"]
        df2 = sdata["transformed_adipocytes_circles"]
        res_symdiff = gpd.overlay(df1, df2, how='difference')
        sdata["new_cell_boundaries_without_adipocytes"] = res_symdiff
        sdata["new_cell_boundaries_without_adipocytes"].index = [f"{self.sample_key}_{i}" for i in range(len(sdata["new_cell_boundaries_without_adipocytes"]))]
        sdata["new_cell_boundaries_with_adipocytes"] = pd.concat([sdata["new_cell_boundaries_without_adipocytes"], sdata["transformed_adipocytes_circles"]])
        sdata["new_cell_boundaries_with_full_adipocytes"] = pd.concat([sdata["new_cell_boundaries_without_adipocytes"], sdata["transformed_adipocytes"]])
        
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
    
    def add_adipo_table(
        self,
        zarr_path=None,
        sdata=None,
        save=False,
        postfix = 'with_adipo_table'
        ):

        """
        Add adipocytes to adata.obs of SpatialData table.
        Loads the existing `new_cell_boundaries_with_adipocytes`
        rebuilds the cell indices (adipocytes contains `_A` in IDs)
        builds a zero-filled AnnData carrying matching metadata,concatenates it with the current table
        update the sdata['table'] with updated meta
        Optionally writes the updated zarr when `save` is True.
        """

    
        # Updating cell index and add sample id along with actual cell id for non adipocytes
        # For adipocytes just use the adipocyte id with sample id appended
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)
        gdf = sdata["new_cell_boundaries_with_adipocytes"]
        fat_indexes = gdf.index[gdf.index.str.contains("_A")].tolist()
        gdf_no_fat = gdf[~gdf.index.isin(fat_indexes)]
        
        # new_cell_boundaries no longer aliggned with the table as in the add adipocytes
        # we thrown away some cells that overlapped with adipocytes
        #cell_gdf = sdata["new_cell_boundaries"] @ Emily code updated
        
        cell_indexes = list(gdf_no_fat.cell_id.astype(str) + f'_{self.sample_key}')
        all_indexes = cell_indexes + fat_indexes
    
        sdata["new_cell_boundaries_with_adipocytes"].index = all_indexes
        sdata["new_cell_boundaries_with_full_adipocytes"].index = all_indexes
        
        # making dummy anndata for fat cells with same meta data as other cells
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
        # fat_df["region_name"] = None
        
        # Creating a dummy transcript matrix with zeros only for adipocytes
        empty_counts = csr_matrix(
            (len(fat_indexes), len(adata.var_names)),
            dtype=np.float32,
        )
        fat_adata = ad.AnnData(empty_counts)
        
        fat_adata.obs_names = fat_indexes
        fat_adata.var_names = list(adata.var_names)
        
        fat_adata.obs = fat_df
        # taking union of two adata as both have same var names
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
                zarr_path = zarr_path.replace(".zarr", f"_{postfix}.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + f"_{postfix}.zarr"), overwrite=True)
        else:
            return sdata
    
    def distance_from_bone(self, zarr_path=None, sdata=None, save=False):
    
    
        sdata, zarr_path = self._resolve_sdata_inputs(zarr_path=zarr_path, sdata=sdata)
            
        sdata["new_cell_boundaries_with_full_adipocytes"].index = sdata["new_cell_boundaries_with_adipocytes"].index.tolist()
        
        bones = sdata["transformed_bone"]
        cells = sdata["new_cell_boundaries"]

        # Getting cells with nearest bone and distance to that bone
        cells_w_bones = gpd.sjoin_nearest(
            cells,
            bones,
            distance_col="distances",
            how="left",
            max_distance=self.CELL_BONE_MAX_DISTANCE
            )
        
        # Cells are likey to be equidistant from multiple bones so keeping one nearest bone per cell
        
        # temp = cells_w_bones.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        # temp = temp.rename(columns={"index_right": "nearest_bone", "distances": "distance_to_bone"})
        # old_meta = sdata.tables["table"].obs.copy()
        # updated_meta = old_meta.merge(temp, left_index=True, right_index=True)
        # updated_meta["bone_region"] = "intertrabecular"
        # updated_meta.loc[(updated_meta["distance_to_bone"] > self.PERITRABECULAR_DISTANCE_THRESHOLD[0]) & (updated_meta["distance_to_bone"] < self.PERITRABECULAR_DISTANCE_THRESHOLD[1]), "bone_region"] = "peritrabecular"
        # updated_meta.loc[updated_meta["distance_to_bone"] < self.ENDOSTEAL_DISTANCE_THRESHOLD, "bone_region"] = "endosteal"
        # from shapely.ops import unary_union
        # union_poly = unary_union(bones.geometry)
        # within_bone = sdata["new_cell_boundaries"].within(union_poly)
        # cell_ids = list(within_bone[within_bone].index)
        # updated_meta.loc[cell_ids, "bone_region"] = "within_bone"

        df_dropped = cells_w_bones.reset_index().drop_duplicates(subset="index", keep='last').set_index("index")
        df_dropped = df_dropped.rename(columns={"index_right": "nearest_bone", "distances": "distance_to_bone"})

        df_dropped = df_dropped[["nearest_bone", "distance_to_bone"]]
        df = sdata.tables["table"].obs.merge(df_dropped, left_index=True, right_index=True)

        sdata.tables["table"].obs["nearest_bone"] = df["nearest_bone"]
        sdata.tables["table"].obs["distance_to_bone"] = df["distance_to_bone"]
        
        sdata.tables["table"].obs["bone_region"] = "intertrabecular"
        sdata.tables["table"].obs.loc[(sdata.tables["table"].obs["distance_to_bone"] > self.PERITRABECULAR_DISTANCE_THRESHOLD[0]) & (sdata.tables["table"].obs["distance_to_bone"] < self.PERITRABECULAR_DISTANCE_THRESHOLD[1]), "bone_region"] = "peritrabecular"
        sdata.tables["table"].obs.loc[sdata.tables["table"].obs["distance_to_bone"] < self.ENDOSTEAL_DISTANCE_THRESHOLD, "bone_region"] = "endosteal"

        #find the cells which are inside a bone polygon

        within_bone = []

        for i in list(bones.index):
            poly = bones.geometry[i]
            outcome = sdata["new_cell_boundaries"].within(poly)
            sub = outcome[outcome == True]
            cell_ids = list(sub.index)
            within_bone.extend(cell_ids)

        #for each cell where distance to bone is 0 call it endosteal
        sdata.tables["table"].obs.loc[within_bone, "bone_region"] = "within_bone"
        if save == True:
            if zarr_path != None:
                zarr_path = zarr_path.replace(".zarr", "_with_bone_distance.zarr")
                sdata.write(zarr_path, overwrite=True)
            else:
                sdata.write(os.path.join(self.save_folder, self.sample_key + "_with_bone_distance.zarr"), overwrite=True)
        else:
            return sdata
            
        

if __name__ == '__main__':
    import subprocess

    cfg = STConfig()
    customizer = SDataCustomizer(cfg, sample_key='10693_R2')
    
    # Instantiating sdata object from xenium run raw data and save it as zar file
    #sdata = customizer.initialise_sdata(save=True)
    
    # Adding H&E Image to sdata using manually obtained landmark matrix from xenium explorer
    # bash_command = customizer.add_he_image()
    # process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()
    # __import__('subprocess').run(bash_command, shell=True, check=True)

    # sdata = customizer.initialise_sdata_from_zarr()
    postfix = 'megs'
    # sdata = customizer.add_megakaryocytes(save=False, postfix=postfix)
    # print(sdata['transformed_megs'].shape)
    # sdata = customizer.add_megs_to_segmentation(
    #     sdata=sdata,
    #     save=True,
    #     postfix=postfix
    #    )
    zarr_path = f'/well/rittscher/users/qwi813/xenium_paper/outputs/10693_R2_{postfix}.zarr'
    #sdata = customizer.initialise_sdata_from_zarr(zarr_path)
    # del sdata.images['07_24_55']
    postfix = 'megs_aggregated'
    #sdata = customizer.aggregate_with_sopa(sdata=sdata, save=True, postfix=postfix)
    # import pdb; pdb.set_trace()

    # zarr_path = f'/well/rittscher/users/qwi813/xenium_paper/outputs/10693_R2_{postfix}.zarr'
    # sdata = customizer.initialise_sdata_from_zarr(zarr_path)
    # sdata = customizer.add_cif_scores(sdata=sdata, save=False, window_size = 512)
    # sdata = customizer.add_meg_phenotypes(sdata=sdata, save=False)
    # sdata = customizer.add_negative_annotations(sdata=sdata, save=False)
    # sdata = customizer.add_positive_annotations(sdata=sdata, save=False)
    postfix = 'with_meta'
    #sdata = customizer.add_metadata(sdata=sdata, save=True, postfix=postfix)
    zarr_path = f'/well/rittscher/users/qwi813/xenium_paper/outputs/10693_R2_{postfix}.zarr'
    sdata = customizer.initialise_sdata_from_zarr(zarr_path)
    sdata = customizer.add_cell_annotations(sdata=sdata, save=False)
    sdata = customizer.change_cell_annotations(sdata=sdata, save=False)
    print('adding bone segmentation')
    sdata = customizer.add_bone_segmentation(sdata=sdata, save=False)
    print('adding adipocytes')
    sdata = customizer.add_adipocytes(sdata=sdata, save=False)
    #postfix = 'with_adipocytes'
    print('adding adipo table')
    sdata = customizer.add_adipo_table(sdata=sdata, save=False, postfix=postfix)
    #zarr_path = f'/well/rittscher/users/qwi813/xenium_paper/outputs/10693_R2_{postfix}.zarr'
    #sdata = customizer.initialise_sdata_from_zarr(zarr_path)
    sdata = customizer.distance_from_bone(sdata=sdata, save=True)
    
    

