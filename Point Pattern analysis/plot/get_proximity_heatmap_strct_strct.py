import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import pdb

#-------------------------------------

want_row_cluster = True
want_col_cluster = True

xl_path = 'path/to/data' #From get_metadata
pvalue_path = 'path/to/permutation_data'
save_path = 'path/to/output'

structure = ['Bone','Fat','Arteriole','Sinusoid']
mpn_grp = ['Normal','ET','PV','MF']

if not os.path.exists(save_path):
    os.makedirs(save_path)

for m in range(len(mpn_grp)):
    xl_path_mpn = '{}/{}.xlsx'.format(xl_path,mpn_grp[m])
    df_temp_strct = []
    df_temp_med = []
    p_ = []
    to_array = []
    for s in range(len(structure)):
        pval_path = '{}/{}/{}/permut.xlsx'.format(pvalue_path,structure[s],mpn_grp[m])
        df_temp_s = pd.read_excel(xl_path_mpn,structure[s])
        p_s = pd.read_excel(pval_path)
        if s == 0:
            p_ = p_s.copy()
        else:
            p_ = np.concatenate((p_,p_s))
        to_array_s = np.full((len(df_temp_s),1),structure[s])
        to_array_s = to_array_s.reshape((len(to_array_s),1))
        if s == 0:
            df_temp_strct = np.array(df_temp_s['structure']).reshape(len(df_temp_s),1)
            df_temp_med = np.array(1-df_temp_s['median']).reshape(len(df_temp_s),1)
            to_array = to_array_s.copy()
        else:
            df_temp_strct = np.concatenate((df_temp_strct,np.array(df_temp_s['structure']).reshape(len(df_temp_s),1)))
            df_temp_med = np.concatenate((df_temp_med,np.array(1-df_temp_s['median']).reshape(len(df_temp_s),1)))
            to_array = np.concatenate((to_array,to_array_s))
    
    df = np.concatenate((df_temp_strct,df_temp_med,to_array),axis=1)
    df = pd.DataFrame(df)
    df.columns = ['Structure From','Rank','Structure To']
    df['Rank'] = df['Rank'].astype(float)
    p_ = pd.DataFrame(p_)
    p_.columns = ['Structure From','P Value','Structure To']
    pivoted_df = df.pivot(index="Structure From", columns="Structure To", values="Rank")
    plt.figure(figsize=(30,30))
    map_ = sns.clustermap(pivoted_df, cmap="coolwarm",row_cluster=want_row_cluster,col_cluster=want_col_cluster)

    strct_numbers = len(np.unique(np.array(p_['Structure To'])))
    
    p_strct = np.zeros((strct_numbers,strct_numbers))
    p_strct = p_strct.astype(float)
    
    for i in range(len(pivoted_df.columns)):
        strct_to_idx = np.where(p_['Structure To'] == pivoted_df.columns[i])[0]
        for ii in range(len(strct_to_idx)):
            strct_from_idx = np.where(pivoted_df.index == p_['Structure From'][strct_to_idx[ii]])[0][0]
            p_strct[strct_from_idx,i] = p_['P Value'][strct_to_idx[ii]]
            
    if want_row_cluster == True:
        for i, ix in enumerate(map_.dendrogram_col.reordered_ind):
            for j, jx in enumerate(map_.dendrogram_row.reordered_ind):
                if jx != ix:
                    text = map_.ax_heatmap.text(
                            i + 0.5,
                            j + 0.5,
                            "*" if p_strct[jx, ix] < 0.05 else "",
                            ha="center",
                            va="center",
                            color="black",
                                                    )
                else:
                    text = map_.ax_heatmap.text(
                            i + 0.5,
                            j + 0.5,
                            "",
                            ha="center",
                            va="center",
                            color="black",
                                                    )
                text.set_fontsize(20)
    else:
        x_ = map_.ax_heatmap.get_xticks()
        y_ = map_.ax_heatmap.get_yticks()
        
        for i in range(len(x_)):
            for ii in range(len(y_)):
                if ii != i:
                    text = map_.ax_heatmap.text(
                            x_[i],
                            y_[ii],
                            "*" if p_strct[ii, i] < 0.05 else "",
                            ha="center",
                            va="center",
                            color="black",
                                                    )
                else:
                    text = map_.ax_heatmap.text(
                            x_[i],
                            y_[ii],
                            "",
                            ha="center",
                            va="center",
                            color="black",
                                                    )
                text.set_fontsize(20)


    plt.savefig("{}/strct_to_strct_heatmap_{}.png".format(save_path,mpn_grp[m]))
    plt.close()
