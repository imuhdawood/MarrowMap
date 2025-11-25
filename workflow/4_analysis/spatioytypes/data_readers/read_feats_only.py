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
        print(cell_of_interest)

        input_file = cfg.pth_spatiotypes_label
        df = pd.read_csv(input_file, index_col='cell_id')
        ids = df.index.tolist()
        features = features = df.columns.to_list()[0:-1]
        dataDf.loc[ids,features] = df.loc[ids, features]
    print(df.shape, dataDf.loc[ids,:].dropna().shape)
    dataDf.to_csv(cfg.pth_spatiotypes_feat_only)
