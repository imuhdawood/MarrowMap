import os
import pandas as pd
import numpy as np
import statistics
import pdb

xl_path = 'path/to/data'

permut_path = 'path/to/permutation'

save_path = 'path/to/output'

cell_xl = []

for file in os.listdir(xl_path):
    cell_xl.append(file)

numb_cell = len(cell_xl)

header = ['Cell Type From', 'Cell Type To','Rank','P Value']

# header_permut = ['Cell Type From', 'P Value','Cell Type To']

xl_temp = pd.read_excel('{}/{}'.format(xl_path,cell_xl[0]),'All')

cell_type_list = np.array(xl_temp['Cell Type'])

# df_all = np.zeros((len(xl_temp['Cell Type'])*len(xl_temp['Cell Type']),4))
# df_all = df_all.astype(object)

df_n = np.zeros((len(xl_temp['Cell Type'])*len(xl_temp['Cell Type']),4))
df_n = df_n.astype(object)

df_et = np.zeros((len(xl_temp['Cell Type'])*len(xl_temp['Cell Type']),4))
df_et = df_et.astype(object)

df_pv = np.zeros((len(xl_temp['Cell Type'])*len(xl_temp['Cell Type']),4))
df_pv = df_pv.astype(object)

df_mf = np.zeros((len(xl_temp['Cell Type'])*len(xl_temp['Cell Type']),4))
df_mf = df_mf.astype(object)

for i in range(len(cell_type_list)):
    # df_all[i*len(cell_type_list):(i+1)*len(cell_type_list),0] = cell_type_list
    # df_all[i*len(cell_type_list):(i+1)*len(cell_type_list),1] = np.full(len(cell_type_list), cell_type_list[i])
    
    df_n[i*len(cell_type_list):(i+1)*len(cell_type_list),0] = cell_type_list
    df_n[i*len(cell_type_list):(i+1)*len(cell_type_list),1] = np.full(len(cell_type_list), cell_type_list[i])
    
    df_et[i*len(cell_type_list):(i+1)*len(cell_type_list),0] = cell_type_list
    df_et[i*len(cell_type_list):(i+1)*len(cell_type_list),1] = np.full(len(cell_type_list), cell_type_list[i])
    
    df_pv[i*len(cell_type_list):(i+1)*len(cell_type_list),0] = cell_type_list
    df_pv[i*len(cell_type_list):(i+1)*len(cell_type_list),1] = np.full(len(cell_type_list), cell_type_list[i])
    
    df_mf[i*len(cell_type_list):(i+1)*len(cell_type_list),0] = cell_type_list
    df_mf[i*len(cell_type_list):(i+1)*len(cell_type_list),1] = np.full(len(cell_type_list), cell_type_list[i])

for i in range(len(cell_xl)):
    cell_type = cell_xl[i].split('.')[0]
    # xl_all = pd.read_excel('{}/{}'.format(xl_path,cell_xl[i]),'All')
    xl_n = pd.read_excel('{}/{}'.format(xl_path,cell_xl[i]),'Normal')
    xl_et = pd.read_excel('{}/{}'.format(xl_path,cell_xl[i]),'ET')
    xl_pv = pd.read_excel('{}/{}'.format(xl_path,cell_xl[i]),'PV')
    xl_mf = pd.read_excel('{}/{}'.format(xl_path,cell_xl[i]),'MF')
    
    if cell_type != 'Granulocyte_mast':
        # xl_all_permut = pd.read_excel('{}/{}/All/permut.xlsx'.format(permut_path,cell_type))
        xl_n_permut = pd.read_excel('{}/{}/Normal/permut.xlsx'.format(permut_path,cell_type))
        xl_et_permut = pd.read_excel('{}/{}/ET/permut.xlsx'.format(permut_path,cell_type))
        xl_pv_permut = pd.read_excel('{}/{}/PV/permut.xlsx'.format(permut_path,cell_type))
        xl_mf_permut = pd.read_excel('{}/{}/MF/permut.xlsx'.format(permut_path,cell_type))
    else:
        # xl_all_permut = pd.read_excel('{}/Granulocyte/mast/All/permut.xlsx'.format(permut_path))
        xl_n_permut = pd.read_excel('{}/Granulocyte_mast/Normal/permut.xlsx'.format(permut_path))
        xl_et_permut = pd.read_excel('{}/Granulocyte_mast/ET/permut.xlsx'.format(permut_path))
        xl_pv_permut = pd.read_excel('{}/Granulocyte_mast/PV/permut.xlsx'.format(permut_path))
        xl_mf_permut = pd

    for ii in range(len(xl_n['Cell Type'])):
        if cell_type == 'Granulocyte_mast':
            cell_type = 'Granulocyte/mast'
        cell_from = xl_n['Cell Type'][ii]
        cell_idx = np.where((np.array(df_n[:,1]) == cell_type) & (np.array(df_n[:,0]) == cell_from))[0][0]
        df_n[cell_idx,2] = statistics.median(np.array(xl_n.iloc[ii][1:]))
        permut_cell_idx = np.where(np.array(xl_n_permut.iloc[:,1]) == cell_from)[0][0]
        df_n[cell_idx,3] = xl_n_permut.iloc[:,2][permut_cell_idx]
    
    for ii in range(len(xl_et['Cell Type'])):
        if cell_type == 'Granulocyte_mast':
            cell_type = 'Granulocyte/mast'
        cell_from = xl_n['Cell Type'][ii]
        cell_from = xl_et['Cell Type'][ii]
        cell_idx = np.where((np.array(df_et[:,1]) == cell_type) & (np.array(df_et[:,0]) == cell_from))[0][0]
        df_et[cell_idx,2] = statistics.median(np.array(xl_et.iloc[ii][1:]))
        permut_cell_idx = np.where(np.array(xl_et_permut.iloc[:,1]) == cell_from)[0][0]
        df_et[cell_idx,3] = xl_et_permut.iloc[:,2][permut_cell_idx]
    
    for ii in range(len(xl_pv['Cell Type'])):
        if cell_type == 'Granulocyte_mast':
            cell_type = 'Granulocyte/mast'
        cell_from = xl_pv['Cell Type'][ii]
        cell_idx = np.where((np.array(df_pv[:,1]) == cell_type) & (np.array(df_pv[:,0]) == cell_from))[0][0]
        df_pv[cell_idx,2] = statistics.median(np.array(xl_pv.iloc[ii][1:]))
        permut_cell_idx = np.where(np.array(xl_pv_permut.iloc[:,1]) == cell_from)[0][0]
        df_pv[cell_idx,3] = xl_pv_permut.iloc[:,2][permut_cell_idx]
    
    for ii in range(len(xl_mf['Cell Type'])):
        if cell_type == 'Granulocyte_mast':
            cell_type = 'Granulocyte/mast'
        cell_from = xl_mf['Cell Type'][ii]
        cell_idx = np.where((np.array(df_mf[:,1]) == cell_type) & (np.array(df_mf[:,0]) == cell_from))[0][0]
        df_mf[cell_idx,2] = statistics.median(np.array(xl_mf.iloc[ii][1:]))
        permut_cell_idx = np.where(np.array(xl_mf_permut.iloc[:,1]) == cell_from)[0][0]
        df_mf[cell_idx,3] = xl_mf_permut.iloc[:,2][permut_cell_idx]
    
    
df_n = pd.DataFrame(df_n)
df_n.columns = header
df_n['Cell Type From'] = df_n['Cell Type From'].astype(str)
df_n['Rank'] = df_n['Rank'].astype(float)
df_n['P Value'] = df_n['P Value'].astype(float)
df_n['Cell Type To'] = df_n['Cell Type To'].astype(str)


df_et = pd.DataFrame(df_et)
df_et.columns = header
df_et['Cell Type From'] = df_et['Cell Type From'].astype(str)
df_et['Rank'] = df_et['Rank'].astype(float)
df_et['P Value'] = df_et['P Value'].astype(float)
df_et['Cell Type To'] = df_et['Cell Type To'].astype(str)


df_pv = pd.DataFrame(df_pv)
df_pv.columns = header
df_pv['Cell Type From'] = df_pv['Cell Type From'].astype(str)
df_pv['Rank'] = df_pv['Rank'].astype(float)
df_pv['P Value'] = df_pv['P Value'].astype(float)
df_pv['Cell Type To'] = df_pv['Cell Type To'].astype(str)

df_mf = pd.DataFrame(df_mf)
df_mf.columns = header
df_mf['Cell Type From'] = df_mf['Cell Type From'].astype(str)
df_mf['Rank'] = df_mf['Rank'].astype(float)
df_mf['Cell Type To'] = df_mf['Cell Type To'].astype(str)


writer = pd.ExcelWriter('{}/cell_cell_prox_gathered.xlsx'.format(save_path), engine="xlsxwriter")
df_n.to_excel(writer, sheet_name="Normal", index=False)
df_et.to_excel(writer, sheet_name="ET", index=False)
df_pv.to_excel(writer, sheet_name="PV", index=False)
df_mf.to_excel(writer, sheet_name="MF", index=False)

writer.close()
