import os
import pandas as pd
import numpy as np

xl_path = 'path/to/ppm_result/structure_cellneighbor'
savepath = 'path/to/structure_cellneighbor'

if not os.path.exists(savepath):
    os.makedirs(savepath)

strct_grp = ['Arteriole','Bone','Fat','Sinusoid']

correspond_mpn = pd.read_excel('Correspond_MPN.xlsx')

mpn_list = np.unique(correspond_mpn['MPN'])

mpn_list = mpn_list[mpn_list != 'PrePMF'] 

cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']

for g in range(0,len(strct_grp)):
    struct_name = strct_grp[g]
    
    xl_list = []
    xl_ids = []

    sub_xl_path = '{}/{}'.format(xl_path,struct_name)

    for file in os.listdir(sub_xl_path):
        xl_list.append(file)
        xl_ids.append(file.split('.csv')[0])
    
    writer_ = pd.ExcelWriter('{}/{}.xlsx'.format(savepath,struct_name), engine="xlsxwriter")
    
    for mp in range(0,len(mpn_list)):
        header = []
        cnt = 0
        list(set(xl_ids) & set(correspond_mpn['ID']))
        col1 = [b for b, c in zip(correspond_mpn['ID'], correspond_mpn['MPN']) if b in list(set(xl_ids) & set(correspond_mpn['ID']))]
        col2 = [c for b, c in zip(correspond_mpn['ID'], correspond_mpn['MPN']) if b in list(set(xl_ids) & set(correspond_mpn['ID']))]

        per_mpn = np.column_stack((col1, col2))
        numb_mpn = len(np.where(np.array(per_mpn[:,1]) == mpn_list[mp])[0])
        x_data_cluster = np.zeros((len(cluster_type_list),numb_mpn+1))
        x_data_cluster = x_data_cluster.astype(object)
        x_data_cluster[:,0] = cluster_type_list
        
        for i in range(len(xl_list)):
            dataname = xl_list[i].split('.csv')[0]
            mpn_type = per_mpn[:,1][np.where(np.array(per_mpn[:,0] == dataname))[0][0]]
            if mpn_type == mpn_list[mp]:
                header.append('{} {}'.format(dataname,mpn_list[mp]))
        
        for i in range(len(xl_list)):
            dataname = xl_list[i].split('.csv')[0]
            mpn_type = per_mpn[:,1][np.where(np.array(per_mpn[:,0] == dataname))[0][0]]
            if mpn_type != mpn_list[mp]:
                continue
            xl_whole = pd.read_csv('{}/{}/{}'.format(xl_path,struct_name,xl_list[i]))
            sort_ = np.argsort(xl_whole['Distance to {}'.format(struct_name)])
            
            for ii in range(len(sort_)):
                type_idx = np.where(np.array(cluster_type_list) == str(xl_whole['Cluster Type'][sort_[ii]]))[0][0]
                rank_ = ii/max(sort_)
                x_data_cluster[type_idx,cnt+1] = rank_
            cnt+=1
            
        df_cell = pd.DataFrame(x_data_cluster)
        df_cell = df_cell.rename(columns={df_cell.columns[0]: 'Cluster Type'})
        df_cell.set_index('Cluster Type', inplace=True)
        df_cell.columns = header
        medians_ = df_cell.iloc[:, 0:].median(axis=1)
        df_cell['median'] = medians_
        df_cell.to_excel(writer_, sheet_name='{}'.format(mpn_list[mp]))
            
    writer_.close()
