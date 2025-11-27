import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

main_data = 'path/to/data'

main_savepth = 'proximity_heatmap'

compare_type = ['struct_struct', 'struct_celltype', 'struct_cellneighbor',
                'celltype_celltype','cellneighbor_cellneighbor']

mpn_grp = ['Normal','ET','PV','MF']

want_row_cluster = True
want_col_cluster = True

for comptype in compare_type:
    if comptype == 'struct_struct':
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
        for m in range(len(mpn_grp)):
            prox_art = pd.read_excel('{}/Prox_arteriole.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
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
            
            plt.savefig("{}/{}.png".format(save_folder,mpn_grp[m]))
            plt.close()
        
    elif comptype == 'struct_celltype':
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
        for m in range(len(mpn_grp)):
            prox_art = pd.read_excel('{}/Prox_arteriole.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
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
            
            plt.savefig("{}/{}.png".format(save_folder,mpn_grp[m]))
            plt.close()
    elif comptype == 'struct_cellneighbor':
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
        for m in range(len(mpn_grp)):
            prox_art = pd.read_excel('{}/Prox_arteriole.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
            prox_bone = pd.read_excel('{}/Prox_Bone.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
            prox_fat = pd.read_excel('{}/Prox_Fat.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
            prox_sinu = pd.read_excel('{}/Prox_Sinusoid_filtered.xlsx'.format(xl_path),'{}'.format(mpn_grp[m]))
            
            # pdb.set_trace()
            
            from_array_art = np.array(prox_art['Cluster Type'])
            from_array_art = from_array_art.reshape((len(from_array_art),1))
            from_array_bone= np.array(prox_bone['Cluster Type'])
            from_array_bone = from_array_bone.reshape((len(from_array_bone),1))
            from_array_fat = np.array(prox_fat['Cluster Type'])
            from_array_fat = from_array_fat.reshape((len(from_array_fat),1))
            from_array_sinu = np.array(prox_sinu['Cluster Type'])
            from_array_sinu = from_array_sinu.reshape((len(from_array_sinu),1))
            
            from_array = np.concatenate((from_array_art,from_array_bone,
                                         from_array_fat,from_array_sinu))
            
            numb_cell = np.unique(from_array)
            
            # pdb.set_trace()
            
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
            df.columns = ['Cluster Type', 'Proximity Rank', 'Structure']
            df['Proximity Rank'] = df['Proximity Rank'].astype(float)
            
            p_df = np.concatenate((from_array,p_,to_array),axis=1)
            p_df = pd.DataFrame(p_df)
            p_df.columns = ['Cluster Type','P Value','Structure']
            p_df['P Value'] = p_df['P Value'].astype(float)
            
            pivoted_df = df.pivot(index="Cluster Type", columns="Structure", values="Proximity Rank")
            
            order = ["Bone", "Fat", "Arteriole", "Sinusoid"]
            pivoted_df = pivoted_df[order]
            plt.figure(figsize=(30,30))
            # map_ = sns.clustermap(pivoted_df, cmap="coolwarm",row_cluster=True,col_cluster=True,vmin=0, vmax=1)
            
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
                    strct_from_idx = np.where(pivoted_df.index == p_df['Cluster Type'][strct_to_idx[ii]])[0][0]
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
            
            plt.savefig("{}/{}.png".format(save_folder,mpn_grp[m]))
            plt.close()
            
    elif comptype == 'celltype_celltype':
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
        for m in range(len(mpn_grp)):
            df = pd.read_excel(xl_path,mpn_grp[m])
            df['Rank'] = 1-df['Rank']
            
            pivoted_df = df.pivot(index="Cell Type From", columns="Cell Type To", values="Rank")

            plt.figure(figsize=(30,30))
            map_ = sns.clustermap(pivoted_df, cmap="coolwarm",row_cluster=want_row_cluster,col_cluster=want_col_cluster)

            cell_numbers = len(np.unique(np.array(df['Cell Type To'])))

            p_cell = np.zeros((cell_numbers,cell_numbers))
            p_cell = p_cell.astype(float)

            for i in range(len(pivoted_df.columns)):
                # pdb.set_trace()
                if pivoted_df.columns[i] != 'Granulocyte/mast':
                    cell_to_idx = np.where(df['Cell Type To'] == pivoted_df.columns[i])[0]
                else:
                    cell_to_idx = np.where(df['Cell Type To'] == 'Granulocyte/mast')[0]
                    
                for ii in range(len(cell_to_idx)):
                    cell_from_idx = np.where(pivoted_df.index == df['Cell Type From'][cell_to_idx[ii]])[0][0]
                    p_cell[cell_from_idx,i] = df['P Value'][cell_to_idx[ii]]

            if want_row_cluster == True:
                for i, ix in enumerate(map_.dendrogram_col.reordered_ind):
                    for j, jx in enumerate(map_.dendrogram_row.reordered_ind):
                        if jx != ix:
                            text = map_.ax_heatmap.text(
                                    i + 0.5,
                                    j + 0.5,
                                    "*" if p_cell[jx, ix] < 0.05 else "",
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
                        text = map_.ax_heatmap.text(
                                x_[i],
                                y_[ii],
                                "*" if p_cell[ii, i] < 0.05 else "",
                                ha="center",
                                va="center",
                                color="black",
                                                        )
                        text.set_fontsize(20)

            plt.savefig("{}/{}.png".format(save_folder,mpn_grp[m]))
            plt.close()
    
    elif comptype == 'cellneighbor_cellneighbor':
        
        cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']
        
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
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
                
                plt.savefig("{}/{}.png".format(save_folder,mpn_grp[m]))
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
            
                plt.savefig("{}/{}_nocluster.png".format(save_folder,mpn_grp[m]))
                plt.close()
        
