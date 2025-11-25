import pandas as pd

import sys
sys.path.append('.')

from utils.utils_io import mkdir
from utils.misc import celltypes_feats
from configuration import STConfig

cfg = STConfig()
PATH = f'{cfg.pth_downstream_out}/obj_w_Cell_annotations_RC_22122024_final_with_osteo/BandFeatures'
OUTPATH = f'{cfg.pth_downstream_out}/obj_w_Cell_annotations_RC_22122024_final_with_osteo/BandFeatureCluster'
from tqdm import tqdm

if __name__ == "__main__":
    dataDf  = pd.read_csv(cfg.pth_cell_annotations_final, index_col='Unnamed: 0')
    for cell_of_interest in tqdm(celltypes_feats):
        file_id = 'Granulocyte_mast' if cell_of_interest in ['Granulocyte/mast'] else cell_of_interest 
        print('Processing *******************************************')
        print('********************************************************')
        print(file_id)
        input_file = f'{OUTPATH}/{file_id}/csv/{file_id}.csv'
        df = pd.read_csv(input_file, index_col='cell_id')
        df['ENVIRON'] = f"{cell_of_interest}-"+df['cluster'].astype(str)
        ids = df.index.tolist()
        dataDf.loc[ids,'ENVIRON'] = df.loc[ids, 'ENVIRON']
    print(df.shape, dataDf.loc[ids, 'ENVIRON'].shape)
    dataDf.to_csv(cfg.pth_spatiotypes_label)
