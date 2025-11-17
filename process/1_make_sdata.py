import sdata_initialiser
from sdata_initialiser import SDATA_INITIALISER, SDATA_CUSTOMISER #importing bespoke functions and classes
import os
import subprocess

#defining a folder to save zarr objects in
save_folder = "<name desired folder>"

for i in range(1, 33): #going sample by sample
    
    sdata = SDATA_CUSTOMISER(num=i, save_folder="<name desired folder> (same as above)") #load data
    
    obj = sdata.initialise_sdata(save=True) #initialize zarr object
    
    bashCommand = sdata.add_he_image() #adding in HE image
    print(bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    file_name = os.path.join(save_folder, sdata.sample_key + ".zarr") #assign the filepath for the zarr object to a variable
    
    obj = sdata.add_megakaryocytes(zarr_path=file_name, save=False) #add megakaryocytes to object
    
    sdata.add_megs_to_segmentation(zarr_path=None, sdata=obj, save=True) #add megakaryocytes to segmentation mask
    
    file_name = os.path.join(save_folder, sdata.sample_key + "_with_new_seg.zarr") #assign the filepath for the zarr object with new segmentations (including megs) to a variable
    
    obj = sdata.initialise_sdata_from_zarr(file_name) #initialize new object containing updated segmentations
    
    sdata.aggregate_with_sopa(zarr_path=None, sdata=obj, save=True) #aggregate with sopa (for concept, see: https://gustaveroussy.github.io/sopa/tutorials/comseg/#4-aggregation:~:text=4.-,Aggregation,-%C2%B6)
    
    file_name = os.path.join(save_folder, sdata.sample_key + "_sopa_aggregated.zarr") #assign the filepath for the zarr object with new segmentations (including megs), aggregated with sopa to a variable
    
    obj = sdata.initialise_sdata_from_zarr(file_name) #initialize new object containing updated segmentations + sopa aggregation
    
    obj = sdata.add_512_cif_scores(zarr_path=None, sdata=obj, save=False, cif_path="<insert path containing input CIF scores>", save_transformed_cif = "<define a path to save the transformed CIF scores to>") #add CIF scores for 512 tile size
    
    obj = sdata.add_1024_cif_scores(zarr_path=None, sdata=obj, save=False, cif_path="<insert path containing input CIF scores>", save_transformed_cif = "<define a path to save the transformed CIF scores to>") #add CIF scores for 1024 tile size
    
    obj = sdata.add_negative_annotations(zarr_path=None, sdata=obj, save=False) #add negative annotation information to object
    
    obj = sdata.add_positive_annotations(zarr_path=None, sdata=obj, save=False) #add positive annotation information to object
    
    obj = sdata.add_metadata(zarr_path=None, sdata=obj, save=False) #add metadata
    
    obj = sdata.add_cell_annotations(zarr_path=None, sdata=obj, save=False, anno_path="<add path for cell annotations>", original_cluster_name="SCT_snn_res.0.3") #add cell annotations (this may be done retrospectively once seurat has been done and cell annotations have been prepared.)
    
    obj = sdata.change_cell_annotations(zarr_path=None, sdata=obj, save=False, map_path = "<add path for revised cell annotations>", original_cluster_name="SCT_snn_res.0.3") #change cell annotations (this may be done retrospectively once seurat has been done and cell annotations have been prepared.)
    
    obj = sdata.add_cif_abundance(zarr_path=None, sdata=obj, save=False, tile_path = "<add path for transformed CIF scores for 512 tile size>", tile_size = 512) #add 512 CIF scores to object
    
    obj = sdata.add_cif_abundance(zarr_path=None, sdata=obj, save=False, tile_path = "<add path for transformed CIF scores for 1024 tile size>", tile_size = 1024) #add 1024 CIF scores to object
    
    obj = sdata.add_poly_abundance(zarr_path=None, sdata=obj, save=False)
    
    obj = sdata.add_bone_segmentation(zarr_path=None, sdata=obj, save=False) #insert bone segmentation
    
    obj = sdata.add_adipocytes(zarr_path=None, sdata=obj, save=False) #insert fat segmentations
    
    del obj.images[sdata.he_key]     #delete he image
    
    obj.write(os.path.join(save_folder, sdata.sample_key + "_no_he.zarr"), overwrite=True) #save final processed object into usable zarr and proceed to next steps.
