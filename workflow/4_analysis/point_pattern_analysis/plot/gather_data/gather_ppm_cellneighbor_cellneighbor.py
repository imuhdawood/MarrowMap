import os
import pandas as pd
import numpy as np

xl_path = 'path/to/ppm_result/cellneighbor_cellneighbor' #From 
savepath = 'path/to/cellneighbor_cellneighbor'

if not os.path.exists(savepath):
    os.makedirs(savepath)

correspond_mpn_all = pd.read_excel('Correspond_MPN.xlsx')

cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']

header_MPN = ['ID','MPN']

for r in range(len(cluster_type_list)):
    header_n = ['Cluster Type']
    header_et = ['Cluster Type']
    header_pv = ['Cluster Type']
    header_mf = ['Cluster Type']
    cluster_to_name = cluster_type_list[r]
    
    sub_xl_folder = '{}/{}'.format(xl_path,cluster_to_name)
    sub_xl_list = os.listdir(sub_xl_folder)
    
    correspond_mpn = np.zeros((len(sub_xl_list),2))
    correspond_mpn = correspond_mpn.astype(object)
    for i in range(len(correspond_mpn)):
        correspond_mpn[i,0] = sub_xl_list[i].split('.')[0]
        idx = np.where(np.array(correspond_mpn_all) == sub_xl_list[i].split('.')[0])[0][0]
        correspond_mpn[i,1] = correspond_mpn_all['MPN'][idx]
    
    correspond_mpn = pd.DataFrame(correspond_mpn)
    correspond_mpn.columns = header_MPN
    
    mpn_list = np.unique(correspond_mpn['MPN'])

    for mp in mpn_list:
        numb_mpn = len(np.where(np.array(correspond_mpn['MPN']) == mp)[0])
        if mp =='Normal':
            x_data_cluster_n = np.zeros((len(cluster_type_list),numb_mpn+1))
            x_data_cluster_n = x_data_cluster_n.astype(object)
            x_data_cluster_n[:,0] = cluster_type_list
        elif mp =='ET':
            x_data_cluster_et = np.zeros((len(cluster_type_list),numb_mpn+1))
            x_data_cluster_et = x_data_cluster_et.astype(object)
            x_data_cluster_et[:,0] = cluster_type_list
        elif mp =='PV':
            x_data_cluster_pv = np.zeros((len(cluster_type_list),numb_mpn+1))
            x_data_cluster_pv = x_data_cluster_pv.astype(object)
            x_data_cluster_pv[:,0] = cluster_type_list
        elif mp =='MF':
            x_data_cluster_mf = np.zeros((len(cluster_type_list),numb_mpn+1))
            x_data_cluster_mf = x_data_cluster_mf.astype(object)
            x_data_cluster_mf[:,0] = cluster_type_list
    
    cnt_n = 0
    cnt_et = 0
    cnt_pv = 0
    cnt_mf = 0
    cnt_prmf = 0

    for i in range(len(sub_xl_list)):
        dataname = sub_xl_list[i].split('.xlsx')[0]
        mpn_type = correspond_mpn['MPN'][np.where(np.array(correspond_mpn['ID'] == dataname))[0][0]]
        
        if mpn_type == 'Normal':
            header_n.append('{} {}'.format(dataname,mpn_type))
        elif mpn_type == 'ET':
            header_et.append('{} {}'.format(dataname,mpn_type))
        elif mpn_type == 'PV':
            header_pv.append('{} {}'.format(dataname,mpn_type))
        elif mpn_type == 'MF':
            header_mf.append('{} {}'.format(dataname,mpn_type))
        
        xl_whole = pd.read_excel('{}/{}'.format(sub_xl_folder,sub_xl_list[i]))
        
        sort_ = np.argsort(xl_whole['distance_to'])
        
        for ii in range(len(sort_)):
            type_idx = np.where(np.array(cluster_type_list) == str(xl_whole['name'][sort_[ii]]))[0][0]
            rank_ = ii/max(sort_)
            
            if mpn_type == 'Normal':
                x_data_cluster_n[type_idx,cnt_n+1] = rank_
            elif mpn_type == 'ET':
                x_data_cluster_et[type_idx,cnt_et+1] = rank_
            elif mpn_type == 'PV':
                x_data_cluster_pv[type_idx,cnt_pv+1] = rank_
            elif mpn_type == 'MF':
                x_data_cluster_mf[type_idx,cnt_mf+1] = rank_
        
        if mpn_type == 'Normal':
            cnt_n+=1
        elif mpn_type == 'ET':
            cnt_et+=1
        elif mpn_type == 'PV':
            cnt_pv+=1
        elif mpn_type == 'MF':
            cnt_mf+=1
    
    writer = pd.ExcelWriter('{}/cluster_{}.xlsx'.format(savepath,cluster_to_name), engine="xlsxwriter")
    
    
    df_cluster_n = pd.DataFrame(x_data_cluster_n)
    df_cluster_n.columns = header_n
    medians_n = df_cluster_n.iloc[:, 1:].median(axis=1)
    df_cluster_n['median'] = medians_n
    df_cluster_n.set_index('Cluster Type', inplace=True)
    df_cluster_n.to_excel(writer, sheet_name="Normal")
    
    df_cluster_et = pd.DataFrame(x_data_cluster_et)
    df_cluster_et.columns = header_et
    medians_et = df_cluster_et.iloc[:, 1:].median(axis=1)
    df_cluster_et['median'] = medians_et
    df_cluster_et.set_index('Cluster Type', inplace=True)
    df_cluster_et.to_excel(writer, sheet_name="ET")
    
    df_cluster_pv = pd.DataFrame(x_data_cluster_pv)
    df_cluster_pv.columns = header_pv
    medians_pv = df_cluster_pv.iloc[:, 1:].median(axis=1)
    df_cluster_pv['median'] = medians_pv
    df_cluster_pv.set_index('Cluster Type', inplace=True)
    df_cluster_pv.to_excel(writer, sheet_name="PV")
    
    df_cluster_mf = pd.DataFrame(x_data_cluster_mf)
    df_cluster_mf.columns = header_mf
    medians_mf = df_cluster_mf.iloc[:, 1:].median(axis=1)
    df_cluster_mf['median'] = medians_mf
    df_cluster_mf.set_index('Cluster Type', inplace=True)
    df_cluster_mf.to_excel(writer, sheet_name="MF")
    writer.close()
