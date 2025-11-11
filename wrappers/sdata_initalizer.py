"""
Instantiate and save a spatialdata.SpatialData object from a Xenium experiment.

This module provides SDataInitalizer, which selects a sample (by sample_key or
row_idx) from a metadata CSV (STConfig) and constructs a xenium SpatialData
object. Use initialise_sdata(...) to read and optionally write a zarr archive,
or initialise_sdata_from_zarr(...) to load an existing zarr.

Usage (example):
    from configuration import STConfig
    cfg = STConfig()
    init = SDataInitalizer(cfg, sample_key="18669_R1")
    sdata = init.initialise_sdata(save=True)

Notes:
- STConfig should be validated before use (ensure paths exist / output dirs
  are writable).
- The script is intended for interactive and scripted workflows; see the
  __main__ block for a simple runnable example.
"""

import os
import sys
from pathlib import Path
import pandas as pd
import spatialdata as sd
from spatialdata_io import xenium

import sys
sys.path.append('/well/rittscher/users/qwi813/xenium_paper')
from configuration import STConfig

class SDataInitalizer:
    def __init__(self, config,
                row_idx =0,
                sample_key = None,
                key_col = 'sample_key'
                 ):
        # import pdb;pdb.set_trace()
        self.keys_csv = pd.read_csv(config.pth_meta_csv)
        self.save_folder = config.pth_out_zar
        self.annotation_path = config.pth_annotation
        self.cif_path = config.pth_cif_score
        self.neg_annotation_save = config.pth_neg_annotation_save

        sample_ids = self.keys_csv[key_col].tolist()
        if sample_key and sample_key not in sample_ids:
            raise ValueError(f"sample_key '{sample_key}' not found in {config.pth_meta_csv}")
        else:
            self.sample_key = sample_key
        
        # If row_idx is passed
        if not self.sample_key:
            if row_idx < 0 or row_idx > len(sample_ids):
                raise IndexError(f"num index {row_idx} out of range (0..{len(self.keys_csv)-1})")
            else:
                self.sample_key = sample_ids[row_idx]
        
        # populate attributes from the selected row
        meta_dict = self.keys_csv[self.keys_csv[key_col] == sample_key].iloc[0].to_dict()
        self.base_path = meta_dict["base_path"]
        self.exp_name = meta_dict["xen_exp"]
        self.he_key = meta_dict["xenium file"]
        self.he_image_path = os.path.join(config.pth_he_img, f"{self.he_key}.ome.tif")
        self.meg_path = meta_dict["meg_path"]
        self.meg_pheno_path = meta_dict["meg_pheno_path"]
        self.mut_status = meta_dict["Mutation_status"]
        self.run = meta_dict["Run"]
        self.study_id = meta_dict["Study_ID"]
        self.diagnosis = meta_dict["Diagnosis"]
        self.diagnosis2 = meta_dict["Diagnosis_2"]
        self.broad_diagnosis = meta_dict["Broad_diagnosis"]
        self.tissue_poly = meta_dict["tissue_poly"]
        
    def initialise_sdata(self, save=True):
        sdata = xenium(os.path.join(self.base_path, self.exp_name))
        # Sdata could be customized here deleting the morphology_mip 
        del sdata.images["morphology_mip"]
        
        if save == True:
            sdata.write(os.path.join(self.save_folder, self.sample_key + ".zarr"), overwrite=True)
        return sdata
    
    def initialise_sdata_from_zarr(self, zarr_path=None):
        if not zarr_path:
            zarr_path = os.path.join(self.save_folder, self.sample_key + ".zarr")
        sdata = sd.read_zarr(zarr_path)
        return sdata

if __name__ == '__main__':
    # pass as dictionary if you want to modify or edit in STConfig
    cfg = STConfig()
    init = SDataInitalizer(cfg, sample_key='10693_R2')
    sdata = init.initialise_sdata(save=True)
    sdata = init.initialise_sdata_from_zarr()