import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.metrics import roc_auc_score
from helphers import mkdirs,  poly_to_area, generate_its_stats
from data_utils import process_intertrabecular_data, create_bags_from_dataframe

import sys
sys.path.append('.')
from utils.misc import environ_feats
from data import MyDataset,DataLoader
from model import Net
from model_wraper import train,test
from configuration import STConfig

cfg = STConfig()

MIL_OUTPUT = f'{cfg.pth_out_zar}/MIL'
ITS_POLY_FILE = cfg.pth_its_poly
META_FILE = cfg.pth_meta_csv
INPUT_FILE = cfg.pth_spatiotypes_feat_label

N_EPOCHS = 200#100#0
OUT_DIR = f'{MIL_OUTPUT}'
ROI_STATS_FILE =  f'{OUT_DIR}/Stats.csv' # ITS Level Spatiotypes abundnace
CSV_DIR = f'{OUT_DIR}/csv/'

mkdirs(OUT_DIR)
mkdirs(CSV_DIR)

skey,dkey, rkey,end_point = 'sample_key', 'diagnosis2', 'it_regions', 'Normal_vs_MF'
ENVIRON = 'ENVIRON'
CELL_TYPE = 'obj.anno_5_w_megs'
CELL_TYPE = 'obj.anno_5_w_megs_w_stromal'

mDf = pd.read_csv(META_FILE,index_col=skey)
mDf.rename(columns={'Diagnosis_2':dkey},inplace=True)
cases = mDf.index.tolist()
#cases = ['10693_R1','18606_R3']

itrDf = process_intertrabecular_data(
    ITS_POLY_FILE,
    INPUT_FILE,
    ROI_STATS_FILE,
    rkey,
    skey,
    dkey,
    ENVIRON,
    CELL_TYPE,
    poly_to_area,
    environ_selected_feats=None
    )
print()
dataDf = itrDf[~itrDf[rkey].isin(['non_intertrabecular'])]
dataDf[rkey] = dataDf[rkey].astype('float').astype('int')
dataDf.index = [f'{dataDf.loc[idx,skey]}_R{dataDf.loc[idx,rkey]}' for idx in itrDf.index]

feats_cols = environ_feats
labels_mapping = {'Normal':-1, 'PrePMF':1, 'ET':1, 'MF':1, 'PV':1}
dataDf[end_point] = dataDf[dkey].replace(labels_mapping)
dataDfc = dataDf.copy()

dataDfR = dataDf[~dataDf[dkey].isin(['Normal','MF'])]
dataDf = dataDf[dataDf[dkey].isin(['Normal','MF'])]
aucs=[]
# Create bags from dataframes
data_nmf = create_bags_from_dataframe(dataDf, feats_cols, labels_col=end_point, bag_id=skey)
data_rest_mpn = create_bags_from_dataframe(dataDfR, feats_cols, labels_col=end_point, bag_id=skey)

group_mapping = {
    "18612_R1": "18612_R5",  # Group for R1 and R4
    "18612_R4": "18612_R5",
    "18612_R2": "18612_R6",  # Group for R2 and R3
    "18612_R3": "18612_R6"
}


data_nmf['group_id'] = np.array([group_mapping[idx] if idx in group_mapping else idx for idx in data_nmf['samples_id']])

# Create a DataFrame with group-level info
group_labels = pd.DataFrame({
    'group_id': data_nmf['group_id'],
    'labels': data_nmf['labels']
}).drop_duplicates()

# Separate the group IDs by label
pos_groups = group_labels[group_labels['labels'] == 1]['group_id'].unique()
neg_groups = group_labels[group_labels['labels'] == -1]['group_id'].unique()

n_pos = len(pos_groups)
n_neg = len(neg_groups)

# Number of bootstrap iterations
N_BOOTSTRAPS = 20# 20  # Adjust as appropriate

wsi_df = pd.DataFrame()
itr_df = pd.DataFrame()

WW = np.array([])
aucs = []

rng = np.random.RandomState()

for boot_idx in range(N_BOOTSTRAPS):

    train_groups = np.concatenate((
        rng.choice(pos_groups, size=n_pos, replace=True),
        rng.choice(neg_groups, size=n_neg, replace=True)
    ))
    test_groups = list(set(group_labels['group_id']) - set(train_groups))
    if dataDf[dataDf[skey].isin(test_groups)][dkey].nunique()<=1:
        aucs.append(np.nan)
        continue # drop this bootstrap run

    print(
        f"Bootstrap {boot_idx} | #Train Groups: {len(train_groups)} | #Test Groups: {len(test_groups)}",
        f" | #Train Set: {len(set(train_groups))} | #Test Set: {len(set(test_groups))}"
         )

    train_idx = np.array([np.where(np.isin(data_nmf['group_id'], group))[0][0] for group in train_groups])
    test_idx = np.where(np.isin(data_nmf['group_id'], test_groups))[0]

    print(f"Train idx: {len(set(train_idx))} | Test idx: {len(set(test_idx))}")
    bags_tr = data_nmf['bags'][train_idx]
    y_tr = data_nmf['labels'][train_idx]

    # Split into positive and negative
    pos_bags = bags_tr[y_tr == 1]
    neg_bags = bags_tr[y_tr == -1]

    # Create DataLoaders
    pos_loader = DataLoader(MyDataset(pos_bags), batch_size=1)
    neg_loader = DataLoader(MyDataset(neg_bags), batch_size=1)

    # Model Initialization and training
    model = Net(len(feats_cols)).cuda()
    optimizer = optim.Adam(model.parameters(), lr=0.01, weight_decay=1e-2)  # or your best LR best decay 1e-2
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    train(model, pos_loader, neg_loader, device, optimizer, epochs=N_EPOCHS)

    # Collect model weights for inspection
    W = model.out.weight[0].detach().cpu().numpy()
    if WW.shape[0] == 0:
        WW = W
    else:
        WW = np.vstack((WW, W))

    # Cohort Level prediction both on MPN and rest cases
    wsi_pred_dict = {}
    itr_pred_dict = {}
    for dataset_name, data_dict in [('MPN', data_nmf), ('Rest_MPN', data_rest_mpn)]:
        if dataset_name == 'MPN':
            test_indices = test_idx  # OOB portion of MPN
        else:
            # For 'Rest_MPN', decide if you want to test on the entire dataset or something else
            test_indices = None
        test(
            model,
            data_dict,
            wsi_pred_dict,  # WSI-level
            itr_pred_dict,  # ITS-level
            fidx=boot_idx,
            test_indices = test_indices
        )

    sDf = pd.DataFrame.from_dict(wsi_pred_dict).T
    iDf = generate_its_stats(itr_pred_dict, feats_cols)

    fold_auc = roc_auc_score(sDf["labels"], sDf["scores"])
    aucs.append(fold_auc)

    print(f"Bootstrap {boot_idx} AUC = {fold_auc}")

    # Combine predictions
    if wsi_df.empty:
        wsi_df = sDf
        itr_df = iDf
    else:
        wsi_df = pd.concat([wsi_df, sDf])
        itr_df = pd.concat([itr_df, iDf])


# Aggregated Score across all ensembles models
wsi_df['fold_idx'] = wsi_df['fold_idx'].astype(int)
itr_df['fold_idx'] = itr_df['fold_idx'].astype(int)

wsi_df.to_csv(f'{CSV_DIR}/bootstraps_sample_pred.csv')
itr_df.to_csv(f'{CSV_DIR}/bootstraps_itr_pred.csv')

wsi_df = (
    wsi_df.groupby(wsi_df.index)
          .agg({'scores': 'median', 'labels': 'first', 'fold_idx': 'first'})
          .reset_index()
)
print(np.nanmean(aucs))


# Visualization informative features bases on effect
plt.figure(figsize=(20, 8))
itr_df[feats_cols].boxplot(); plt.xticks(rotation=45, ha = 'right'); plt.tight_layout()
plt.savefig(f'{CSV_DIR}/effects_plot.png', dpi=300)

# Compute combined AUC
c_auc = roc_auc_score(wsi_df["labels"], wsi_df["scores"])
print("Bootstrap AUCs:", aucs)
print("Bootstrap AUCs Mean:", np.nanmean(aucs))
print("Combined AUC:", c_auc)

# Save outputs
wsi_df = (
    dataDfc.loc[:, [skey, dkey]]
           .drop_duplicates()
           .reset_index(drop=True)
           .set_index(skey)
           .join(wsi_df.set_index('index'))
)
wsi_df.to_csv(f'{CSV_DIR}/sample_pred.csv')
itr_df.to_csv(f'{CSV_DIR}/itr_pred.csv')
pd.DataFrame(aucs).to_csv(f'{CSV_DIR}/aucs.csv')



