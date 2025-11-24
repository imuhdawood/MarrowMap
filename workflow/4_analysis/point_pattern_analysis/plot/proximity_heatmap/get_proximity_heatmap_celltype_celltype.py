import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

#-------------------------------------

want_row_cluster = True
want_col_cluster = True

xl_path = 'path/to/data'
pvalue_path = 'path/to/pvalue'
save_path = 'path/to/output'

if not os.path.exists(save_path):
    os.makedirs(save_path)

df_all = pd.read_excel(xl_path,'All')
df_all['Rank'] = 1-df_all['Rank']
df_n = pd.read_excel(xl_path,'Normal')
df_n['Rank'] = 1-df_n['Rank']
df_et = pd.read_excel(xl_path,'ET')
df_et['Rank'] = 1-df_et['Rank']
df_pv = pd.read_excel(xl_path,'PV')
df_pv['Rank'] = 1-df_pv['Rank']
df_mf = pd.read_excel(xl_path,'MF')
df_mf['Rank'] = 1-df_mf['Rank']
# df_prmf = pd.read_excel(xl_path,'PrePMF')
# df_prmf['Rank'] = 1-df_prmf['Rank']


p_all = pd.read_excel(pvalue_path,sheet_name='All')
p_n = pd.read_excel(pvalue_path,sheet_name='Normal')
p_et = pd.read_excel(pvalue_path,sheet_name='ET')
p_pv = pd.read_excel(pvalue_path,sheet_name='PV')
p_mf = pd.read_excel(pvalue_path,sheet_name='MF')
# p_prmf = pd.read_excel(pvalue_path,sheet_name='PrePMF')

# p_all = pd.read_excel(pvalue_path,sheet_name='All',index_col=0)
# p_n = pd.read_excel(pvalue_path,sheet_name='Normal',index_col=0)
# p_et = pd.read_excel(pvalue_path,sheet_name='ET',index_col=0)
# p_pv = pd.read_excel(pvalue_path,sheet_name='PV',index_col=0)
# p_mf = pd.read_excel(pvalue_path,sheet_name='MF',index_col=0)
# p_prmf = pd.read_excel(pvalue_path,sheet_name='PrePMF',index_col=0)

pivoted_df_all = df_all.pivot(index="Cell Type From", columns="Cell Type To", values="Rank")

plt.figure(figsize=(30,30))
map_all = sns.clustermap(pivoted_df_all, cmap="coolwarm",row_cluster=want_row_cluster,col_cluster=want_col_cluster)

cell_numbers = len(np.unique(np.array(p_all['Cell Type To'])))

p_all_cell = np.zeros((cell_numbers,cell_numbers))
p_all_cell = p_all_cell.astype(float)

for i in range(len(pivoted_df_all.columns)):
    # pdb.set_trace()
    if pivoted_df_all.columns[i] != 'Granulocyte/mast':
        cell_to_idx = np.where(p_all['Cell Type To'] == pivoted_df_all.columns[i])[0]
    else:
        # pdb.set_trace()
        cell_to_idx = np.where(p_all['Cell Type To'] == 'Granulocyte/mast')[0]
    for ii in range(len(cell_to_idx)):
        cell_from_idx = np.where(pivoted_df_all.index == p_all['Cell Type From'][cell_to_idx[ii]])[0][0]
        p_all_cell[cell_from_idx,i] = p_all['P Value'][cell_to_idx[ii]]
        

if want_row_cluster == True:
    for i, ix in enumerate(map_all.dendrogram_col.reordered_ind):
        for j, jx in enumerate(map_all.dendrogram_row.reordered_ind):
            # if i == 0:
            #     if j == 10:
            #         pdb.set_trace()
            # if i == 0:
            #     pdb.set_trace()
            # pdb.set_trace()
            if jx != ix:
                text = map_all.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "*" if p_all_cell[jx, ix] < 0.05 else "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            else:
                text = map_all.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            text.set_fontsize(20)
else:
    x_ = map_all.ax_heatmap.get_xticks()
    y_ = map_all.ax_heatmap.get_yticks()
    
    for i in range(len(x_)):
        for ii in range(len(y_)):
            text = map_all.ax_heatmap.text(
                    x_[i],
                    y_[ii],
                    "*" if p_all_cell[ii, i] < 0.05 else "",
                    ha="center",
                    va="center",
                    color="black",
                                            )
            text.set_fontsize(20)


plt.savefig("{}/Cell_heatmap_all.png".format(save_path))
plt.close()

pivoted_df_n = df_n.pivot(index="Cell Type From", columns="Cell Type To", values="Rank")

plt.figure(figsize=(30,30))
map_n = sns.clustermap(pivoted_df_n, cmap="coolwarm",row_cluster=want_row_cluster,col_cluster=want_col_cluster)

cell_numbers = len(np.unique(np.array(p_n['Cell Type To'])))

p_n_cell = np.zeros((cell_numbers,cell_numbers))
p_n_cell = p_n_cell.astype(float)

for i in range(len(pivoted_df_n.columns)):
    # pdb.set_trace()
    if pivoted_df_n.columns[i] != 'Granulocyte/mast':
        cell_to_idx = np.where(p_n['Cell Type To'] == pivoted_df_n.columns[i])[0]
    else:
        cell_to_idx = np.where(p_n['Cell Type To'] == 'Granulocyte/mast')[0]
    for ii in range(len(cell_to_idx)):
        cell_from_idx = np.where(pivoted_df_n.index == p_n['Cell Type From'][cell_to_idx[ii]])[0][0]
        p_n_cell[cell_from_idx,i] = p_n['P Value'][cell_to_idx[ii]]

if want_row_cluster == True:
    for i, ix in enumerate(map_n.dendrogram_col.reordered_ind):
        for j, jx in enumerate(map_n.dendrogram_row.reordered_ind):
            if jx != ix:
                text = map_n.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "*" if p_n_cell[jx, ix] < 0.05 else "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            else:
                text = map_n.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            text.set_fontsize(20)

else:
    x_ = map_n.ax_heatmap.get_xticks()
    y_ = map_n.ax_heatmap.get_yticks()
    
    for i in range(len(x_)):
        for ii in range(len(y_)):
            text = map_n.ax_heatmap.text(
                    x_[i],
                    y_[ii],
                    "*" if p_n_cell[ii, i] < 0.05 else "",
                    ha="center",
                    va="center",
                    color="black",
                                            )
            text.set_fontsize(20)

plt.savefig("{}/Cell_heatmap_n.png".format(save_path))
plt.close()

pivoted_df_et = df_et.pivot(index="Cell Type From", columns="Cell Type To", values="Rank")

plt.figure(figsize=(30,30))

map_et = sns.clustermap(pivoted_df_et, cmap="coolwarm",row_cluster=want_row_cluster,col_cluster=want_col_cluster)

cell_numbers = len(np.unique(np.array(p_et['Cell Type To'])))

p_et_cell = np.zeros((cell_numbers,cell_numbers))
p_et_cell = p_et_cell.astype(float)

for i in range(len(pivoted_df_et.columns)):
    # pdb.set_trace()
    if pivoted_df_et.columns[i] != 'Granulocyte/mast':
        cell_to_idx = np.where(p_et['Cell Type To'] == pivoted_df_et.columns[i])[0]
    else:
        cell_to_idx = np.where(p_et['Cell Type To'] == 'Granulocyte/mast')[0]
    for ii in range(len(cell_to_idx)):
        cell_from_idx = np.where(pivoted_df_et.index == p_et['Cell Type From'][cell_to_idx[ii]])[0][0]
        p_et_cell[cell_from_idx,i] = p_et['P Value'][cell_to_idx[ii]]

if want_row_cluster == True:
    for i, ix in enumerate(map_et.dendrogram_col.reordered_ind):
        for j, jx in enumerate(map_et.dendrogram_row.reordered_ind):
            if jx != ix:
                text = map_et.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "*" if p_et_cell[jx, ix] < 0.05 else "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            else:
                text = map_et.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            text.set_fontsize(20)
else:
    x_ = map_et.ax_heatmap.get_xticks()
    y_ = map_et.ax_heatmap.get_yticks()
    
    for i in range(len(x_)):
        for ii in range(len(y_)):
            text = map_et.ax_heatmap.text(
                    x_[i],
                    y_[ii],
                    "*" if p_et_cell[ii, i] < 0.05 else "",
                    ha="center",
                    va="center",
                    color="black",
                                            )
            text.set_fontsize(20)
            
plt.savefig("{}/Cell_heatmap_et.png".format(save_path))
plt.close()

pivoted_df_pv = df_pv.pivot(index="Cell Type From", columns="Cell Type To", values="Rank")

plt.figure(figsize=(30,30))

map_pv = sns.clustermap(pivoted_df_pv, cmap="coolwarm",row_cluster=want_row_cluster,col_cluster=want_col_cluster)

cell_numbers = len(np.unique(np.array(p_pv['Cell Type To'])))

p_pv_cell = np.zeros((cell_numbers,cell_numbers))
p_pv_cell = p_pv_cell.astype(float)

for i in range(len(pivoted_df_pv.columns)):
    # pdb.set_trace()
    if pivoted_df_pv.columns[i] != 'Granulocyte/mast':
        cell_to_idx = np.where(p_pv['Cell Type To'] == pivoted_df_pv.columns[i])[0]
    else:
        cell_to_idx = np.where(p_pv['Cell Type To'] == 'Granulocyte/mast')[0]
    for ii in range(len(cell_to_idx)):
        cell_from_idx = np.where(pivoted_df_pv.index == p_pv['Cell Type From'][cell_to_idx[ii]])[0][0]
        p_pv_cell[cell_from_idx,i] = p_pv['P Value'][cell_to_idx[ii]]

if want_row_cluster == True:
    for i, ix in enumerate(map_pv.dendrogram_col.reordered_ind):
        for j, jx in enumerate(map_pv.dendrogram_row.reordered_ind):
            if jx != ix:
                text = map_pv.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "*" if p_pv_cell[jx, ix] < 0.05 else "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            else:
                text = map_pv.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            text.set_fontsize(20)

else:
    x_ = map_pv.ax_heatmap.get_xticks()
    y_ = map_pv.ax_heatmap.get_yticks()
    
    for i in range(len(x_)):
        for ii in range(len(y_)):
            text = map_pv.ax_heatmap.text(
                    x_[i],
                    y_[ii],
                    "*" if p_pv_cell[ii, i] < 0.05 else "",
                    ha="center",
                    va="center",
                    color="black",
                                            )
            text.set_fontsize(20)

plt.savefig("{}/Cell_heatmap_pv.png".format(save_path))
plt.close()

pivoted_df_mf = df_mf.pivot(index="Cell Type From", columns="Cell Type To", values="Rank")
plt.figure(figsize=(30,30))

map_mf = sns.clustermap(pivoted_df_mf, cmap="coolwarm",row_cluster=want_row_cluster,col_cluster=want_col_cluster)

cell_numbers = len(np.unique(np.array(p_mf['Cell Type To'])))

p_mf_cell = np.zeros((cell_numbers,cell_numbers))
p_mf_cell = p_mf_cell.astype(float)

for i in range(len(pivoted_df_mf.columns)):
    # pdb.set_trace()
    if pivoted_df_mf.columns[i] != 'Granulocyte/mast':
        cell_to_idx = np.where(p_mf['Cell Type To'] == pivoted_df_mf.columns[i])[0]
    else:
        cell_to_idx = np.where(p_mf['Cell Type To'] == 'Granulocyte/mast')[0]
    for ii in range(len(cell_to_idx)):
        cell_from_idx = np.where(pivoted_df_mf.index == p_mf['Cell Type From'][cell_to_idx[ii]])[0][0]
        p_mf_cell[cell_from_idx,i] = p_mf['P Value'][cell_to_idx[ii]]

if want_row_cluster == True:
    for i, ix in enumerate(map_mf.dendrogram_col.reordered_ind):
        for j, jx in enumerate(map_mf.dendrogram_row.reordered_ind):
            if jx != ix:
                text = map_mf.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "*" if p_mf_cell[jx, ix] < 0.05 else "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            else:
                text = map_mf.ax_heatmap.text(
                        i + 0.5,
                        j + 0.5,
                        "",
                        ha="center",
                        va="center",
                        color="black",
                                                )
            text.set_fontsize(20)

else:
    x_ = map_mf.ax_heatmap.get_xticks()
    y_ = map_mf.ax_heatmap.get_yticks()
    
    for i in range(len(x_)):
        for ii in range(len(y_)):
            text = map_mf.ax_heatmap.text(
                    x_[i],
                    y_[ii],
                    "*" if p_mf_cell[ii, i] < 0.05 else "",
                    ha="center",
                    va="center",
                    color="black",
                                            )
            text.set_fontsize(20)

plt.savefig("{}/Cell_heatmap_mf.png".format(save_path))
plt.close()
