import os
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

xl_path = 'path/to/data'

save_path = 'path/to/output'

mpn_grp = ['Normal','ET','PV','MF']

if not os.path.exists(save_path):
    os.makedirs(save_path)

for m in range(len(mpn_grp)):
    prox_art = pd.read_excel('{}/Prox_arteriole_filtered.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
    prox_bone = pd.read_excel('{}/Prox_Bone.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
    prox_fat = pd.read_excel('{}/Prox_Fat.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
    prox_sinu = pd.read_excel('{}/Prox_Sinusoid_filtered.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
    
    from_array_art = np.array(prox_art['Cell Type'])
    from_array_art = from_array_art.reshape((len(from_array_art),1))
    from_array_bone= np.array(prox_bone['Cell Type'])
    from_array_bone = from_array_bone.reshape((len(from_array_bone),1))
    from_array_fat = np.array(prox_fat['Cell Type'])
    from_array_fat = from_array_fat.reshape((len(from_array_fat),1))
    from_array_sinu = np.array(prox_sinu['Cell Type'])
    from_array_sinu = from_array_sinu.reshape((len(from_array_sinu),1))
    
    from_array = np.concatenate((from_array_art,from_array_bone,
                                 from_array_fat,from_array_sinu))
    
    numb_cell = np.unique(from_array)
    
    prox_array_art = np.array(prox_art['Proximity Rank'])
    prox_array_art = 1-prox_array_art
    prox_array_art = prox_array_art.reshape((len(prox_array_art),1))
    prox_array_bone= np.array(prox_bone['Proximity Rank'])
    prox_array_bone = 1-prox_array_bone
    prox_array_bone = prox_array_bone.reshape((len(prox_array_bone),1))
    prox_array_fat = np.array(prox_fat['Proximity Rank'])
    prox_array_fat = 1-prox_array_fat
    prox_array_fat = prox_array_fat.reshape((len(prox_array_fat),1))
    prox_array_sinu = np.array(prox_sinu['Proximity Rank'])
    prox_array_sinu = 1-prox_array_sinu
    prox_array_sinu = prox_array_sinu.reshape((len(prox_array_sinu),1))
    
    prox_array = np.concatenate((prox_array_art,prox_array_bone,
                                 prox_array_fat,prox_array_sinu))
    
    p_art = np.array(prox_art['P-Value'])
    p_art = p_art.reshape((len(p_art),1))
    p_bone = np.array(prox_bone['P-Value'])
    p_bone = p_bone.reshape((len(p_bone),1))
    p_fat = np.array(prox_fat['P-Value'])
    p_fat = p_fat.reshape((len(p_fat),1))
    p_sinu = np.array(prox_sinu['P-Value'])
    p_sinu = p_sinu.reshape((len(p_sinu),1))
    
    p_ = np.concatenate((p_art,p_bone,
                         p_fat,p_sinu))
    
    to_array_art = np.full((len(from_array_art),1),'Arteriole')
    to_array_bone= np.full((len(from_array_bone),1),'Bone')
    to_array_fat = np.full((len(from_array_fat),1),'Fat')
    to_array_sinu = np.full((len(from_array_sinu),1),'Sinusoid')
    
    to_array = np.concatenate((to_array_art,to_array_bone,
                               to_array_fat,to_array_sinu))
    
    df = np.concatenate((from_array,prox_array,to_array),axis=1)
    df = pd.DataFrame(df)
    df.columns = ['Cell Type', 'Proximity Rank', 'Structure']
    df['Proximity Rank'] = df['Proximity Rank'].astype(float)
    
    p_df = np.concatenate((from_array,p_,to_array),axis=1)
    p_df = pd.DataFrame(p_df)
    p_df.columns = ['Cell Type','P Value','Structure']
    p_df['P Value'] = p_df['P Value'].astype(float)
    
    pivoted_df = df.pivot(index="Cell Type", columns="Structure", values="Proximity Rank")
    
    order = ["Bone", "Fat", "Arteriole", "Sinusoid"]
    pivoted_df = pivoted_df[order]
    plt.figure(figsize=(30,30))
    
    map_ = sns.clustermap(
    pivoted_df[order],       # reorder columns manually
    cmap="coolwarm",
    row_cluster=True,
    col_cluster=False,       # disable column clustering
    vmin=0, vmax=1
)
    
    # build p_strct (unchanged logic, using the reordered pivoted_df.columns)
    strct_numbers = len(np.unique(np.array(p_df['Structure'])))
    p_strct = np.zeros((len(numb_cell), strct_numbers), dtype=float)
    
    for i in range(len(pivoted_df.columns)):
        strct_to_idx = np.where(p_df['Structure'] == pivoted_df.columns[i])[0]
        for ii in range(len(strct_to_idx)):
            strct_from_idx = np.where(pivoted_df.index == p_df['Cell Type'][strct_to_idx[ii]])[0][0]
            p_strct[strct_from_idx, i] = p_df['P Value'][strct_to_idx[ii]]
    
    # Determine column ordering index array robustly:
    if hasattr(map_.dendrogram_col, "reordered_ind") and map_.dendrogram_col.reordered_ind is not None:
        col_order_idx = map_.dendrogram_col.reordered_ind
    else:
        # when col_cluster=False, use natural order 0..ncols-1
        col_order_idx = list(range(len(pivoted_df.columns)))
    
    # Row ordering (we still cluster rows so this should exist)
    row_order_idx = map_.dendrogram_row.reordered_ind
    
    # Annotate stars
    for i, ix in enumerate(col_order_idx):
        for j, jx in enumerate(row_order_idx):
            txt = "*" if p_strct[jx, ix] < 0.05 else ""
            text = map_.ax_heatmap.text(
                i + 0.5,
                j + 0.5,
                txt,
                ha="center",
                va="center",
                color="black"
            )
            text.set_fontsize(20)
    
    plt.savefig("{}/{}.png".format(save_path,mpn_grp[m]))
    plt.close()
