import pandas as pd
from tqdm import tqdm

import sys
sys.path.append('.')
from configuration import STConfig
from wrappers.sdata_initalizer import SDataInitalizer
from wrappers.sdata_customizer import SDataCustomizer
import subprocess

cfg = STConfig()
meta_df = pd.read_csv(cfg.pth_meta_csv)
samples = meta_df['sample_key'].tolist()
print(f'Number of samples {len(samples)}')
print(samples)

for skey in tqdm(samples):
    customizer = SDataCustomizer(config=cfg, sample_key=skey)
    # Instantiating sdata object from xenium run raw data and save it as zar file
    sdata = customizer.initialise_sdata(save=True)

    # Adding H&E Image to sdata using manually obtained landmark matrix from xenium explorer
    print("Adding H&E image")
    bash_command = customizer.add_he_image()
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    __import__('subprocess').run(bash_command, shell=True, check=True)
    sdata = customizer.initialise_sdata_from_zarr()
    print("Adding Megs to sdata")
    sdata = customizer.add_megakaryocytes(save=False)
    print(sdata['transformed_megs'].shape)
    print("Modifying Cell Boundaries based on Meg segmentation")
    sdata = customizer.add_megs_to_segmentation(
        sdata=sdata,
        save=False,
       )
    print('Recomputing cell level expression vector based on updated cell boundaries')
    sdata = customizer.aggregate_with_sopa(sdata=sdata, save=False)
    print('Adding CIF Score')
    sdata = customizer.add_cif_scores(sdata=sdata, save=False, window_size = 512)
    print('Adding Meg Phenotype information')
    sdata = customizer.add_meg_phenotypes(sdata=sdata, save=False)
    print('Adding Negative Annotations')
    sdata = customizer.add_negative_annotations(sdata=sdata, save=False)
    print('Adding Positive Annotations')
    sdata = customizer.add_positive_annotations(sdata=sdata, save=False)
    print('Adding Cell Annotations')
    sdata = customizer.add_cell_annotations(sdata=sdata, save=False)
    print('Updating Cell Annotations')
    sdata = customizer.change_cell_annotations(sdata=sdata, save=False)
    print('Adding bone segmentation')
    sdata = customizer.add_bone_segmentation(sdata=sdata, save=False)
    print('Adding adipocytes')
    sdata = customizer.add_adipocytes(sdata=sdata, save=False)
    print('Adding adipo table')
    sdata = customizer.add_adipo_table(sdata=sdata, save=False)
    print('Computing Cells closer to bone')
    sdata = customizer.distance_from_bone(sdata=sdata, save=True)