import os
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

xl_path = 'path/to/data'

save_path = 'path/to/output'

mpn_grp = ['Normal','ET','PV','MF']

cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']

if not os.path.exists(save_path):
    os.makedirs(save_path)

want_row_cluster = False

for m in range(len(mpn_grp)):
    for i in range(len(cluster_type_list)):
        prox_ = pd.read_excel('{}/Prox_{}.xlsx'.format(xl_path,cluster_type_list[i]),'{}'.format(mpn_grp[m]))
        from_array_temp = np.array(prox_['Cluster Type'])
        from_array_temp = from_array_temp.reshape((len(from_array_temp),1))
        prox_array_temp = 1 - np.array(prox_['Proximity Rank'])
        prox_array_temp = prox_array_temp.reshape((len(prox_array_temp),1))
        p_temp = np.array(prox_['P-Value'])
        p_temp = p_temp.reshape((len(p_temp),1))
        to_array_temp = np.full((len(from_array_temp),1),'{}'.format(cluster_type_list[i]))
        
        if i == 0:
            from_array = from_array_temp
            prox_array = prox_array_temp
            p_ = p_temp
            to_array = to_array_temp
        else:
            from_array = np.concatenate((from_array,from_array_temp))
            prox_array = np.concatenate((prox_array,prox_array_temp))
            p_ = np.concatenate((p_,p_temp))
            to_array = np.concatenate((to_array,to_array_temp))
        
    df = np.concatenate((from_array,prox_array,to_array),axis=1)
    df = pd.DataFrame(df)
    df.columns = ['Cluster Type From', 'Proximity Rank', 'Cluster Type To']
    df['Proximity Rank'] = df['Proximity Rank'].astype(float)
    
    # pdb.set_trace()
    
    p_df = np.concatenate((from_array,p_,to_array),axis=1)
    p_df = pd.DataFrame(p_df)
    p_df.columns = ['Cluster Type From','P Value','Cluster Type To']
    p_df['P Value'] = p_df['P Value'].astype(float)
    pivoted_df = df.pivot(index="Cluster Type From", columns="Cluster Type To", values="Proximity Rank")
    
    plt.figure(figsize=(30,30))
    if want_row_cluster == True:
        map_ = sns.clustermap(pivoted_df, cmap="coolwarm",row_cluster=True,col_cluster=True)
    else:
        map_ = sns.clustermap(pivoted_df, cmap="coolwarm",row_cluster=False,col_cluster=False)
    
    numb_cell = np.unique(from_array)
    
    p_strct = np.zeros((len(numb_cell),len(cluster_type_list)))
    
    p_strct = p_strct.astype(float)
    
    for i in range(len(pivoted_df.columns)):
        # pdb.set_trace()
        cluster_to_idx = np.where(p_df['Cluster Type To'] == pivoted_df.columns[i])[0]
        for ii in range(len(cluster_to_idx)):
            cluster_from_idx = np.where(pivoted_df.index == p_df['Cluster Type From'][cluster_to_idx[ii]])[0][0]
            p_strct[cluster_from_idx,i] = p_df['P Value'][cluster_to_idx[ii]]
    
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
        
        plt.savefig("{}/{}.png".format(save_path,mpn_grp[m]))
        plt.close()
    
    else:
        x_ = map_.ax_heatmap.get_xticks()
        y_ = map_.ax_heatmap.get_yticks()
        
        # pdb.set_trace()
        
        for i in range(len(x_)):
            for ii in range(len(y_)):
                if x_[i] != y_[ii]:
                    text = map_.ax_heatmap.text(
                            x_[i],
                            y_[ii],
                            "*" if p_strct[ii, i] < 0.05 else "",
                            ha="center",
                            va="center",
                            color="black",
                                                    )
                    
                    text.set_fontsize(20)
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
    
        plt.savefig("{}/{}_nocluster.png".format(save_path,mpn_grp[m]))
        plt.close()
