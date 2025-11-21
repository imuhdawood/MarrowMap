import os
import pandas as pd
import numpy as np

xl_path = 'path/to/data'

savepath = 'path/to/output'

if not os.path.exists(savepath):
    os.makedirs(savepath)

cell_folder_list = []

for file in os.listdir(xl_path):
    cell_folder_list.append(file)


with open('celltypes.txt') as p:
    strings_p = p.readlines()
    cell_type_list = [string.split()[0] for string in strings_p]
    p.close()

cell_type_list.remove('CD69')

correspond_mpn = pd.read_excel('Correspond_MPN.xlsx')

mpn_list = np.unique(correspond_mpn['MPN'])

for cc in range(len(cell_folder_list)):
    
    xl_list = []

    for file in os.listdir('{}/{}'.format(xl_path,cell_folder_list[cc])):
        xl_list.append(file)
    
    x_data_cell = np.zeros((len(cell_type_list),len(xl_list)+1))
    x_data_cell = x_data_cell.astype(object)
    x_data_cell[:,0] = cell_type_list
    
    for mp in mpn_list:
        numb_mpn = len(np.where(np.array(correspond_mpn['MPN']) == mp)[0])
        if mp =='Normal':
            x_data_cell_n = np.zeros((len(cell_type_list),numb_mpn+1))
            x_data_cell_n = x_data_cell_n.astype(object)
            x_data_cell_n[:,0] = cell_type_list
        elif mp =='ET':
            x_data_cell_et = np.zeros((len(cell_type_list),numb_mpn+1))
            x_data_cell_et = x_data_cell_et.astype(object)
            x_data_cell_et[:,0] = cell_type_list
        elif mp =='PV':
            x_data_cell_pv = np.zeros((len(cell_type_list),numb_mpn+1))
            x_data_cell_pv = x_data_cell_pv.astype(object)
            x_data_cell_pv[:,0] = cell_type_list
        elif mp =='MF':
            x_data_cell_mf = np.zeros((len(cell_type_list),numb_mpn+1))
            x_data_cell_mf = x_data_cell_mf.astype(object)
            x_data_cell_mf[:,0] = cell_type_list
        elif mp =='PrePMF':
            x_data_cell_prmf = np.zeros((len(cell_type_list),numb_mpn+1))
            x_data_cell_prmf = x_data_cell_prmf.astype(object)
            x_data_cell_prmf[:,0] = cell_type_list
    
    header = ['Cell Type']
    header_n = ['Cell Type']
    header_et = ['Cell Type']
    header_pv = ['Cell Type']
    header_mf = ['Cell Type']
    header_prmf = ['Cell Type']
    
    cnt_n = 0
    cnt_et = 0
    cnt_pv = 0
    cnt_mf = 0
    cnt_prmf = 0
    
    for i in range(len(xl_list)):
        dataname = xl_list[i].split('.xlsx')[0]
        mpn_type = correspond_mpn['MPN'][np.where(np.array(correspond_mpn['ID'] == dataname))[0][0]]
        header.append('{} {}'.format(dataname,mpn_type))
        
        if mpn_type == 'Normal':
            header_n.append('{} {}'.format(dataname,mpn_type))
        elif mpn_type == 'ET':
            header_et.append('{} {}'.format(dataname,mpn_type))
        elif mpn_type == 'PV':
            header_pv.append('{} {}'.format(dataname,mpn_type))
        elif mpn_type == 'MF':
            header_mf.append('{} {}'.format(dataname,mpn_type))
        elif mpn_type == 'PrePMF':
            header_prmf.append('{} {}'.format(dataname,mpn_type))
        
        xl_whole = pd.read_excel('{}/{}/{}'.format(xl_path,cell_folder_list[cc],xl_list[i]))
        
        if len(np.where(np.array(xl_whole['name'] == 'CD69'))[0]) > 0:
            c_idx = np.where(np.array(xl_whole['name'] == 'CD69'))[0][0]
            xl_whole = xl_whole.drop([c_idx])
        
        cell_sort = np.argsort(xl_whole['distance_to_cell'])
        
        for ii in range(len(cell_sort)):
            type_idx = np.where(np.array(cell_type_list) == xl_whole['name'][cell_sort[ii]])[0][0]
            cell_rank = ii/max(cell_sort)
            x_data_cell[type_idx,i+1] = cell_rank
            
            if mpn_type == 'Normal':
                x_data_cell_n[type_idx,cnt_n+1] = cell_rank
            elif mpn_type == 'ET':
                x_data_cell_et[type_idx,cnt_et+1] = cell_rank
            elif mpn_type == 'PV':
                x_data_cell_pv[type_idx,cnt_pv+1] = cell_rank
            elif mpn_type == 'MF':
                x_data_cell_mf[type_idx,cnt_mf+1] = cell_rank
            elif mpn_type == 'PrePMF':
                x_data_cell_prmf[type_idx,cnt_prmf+1] = cell_rank
        
        if mpn_type == 'Normal':
            cnt_n+=1
        elif mpn_type == 'ET':
            cnt_et+=1
        elif mpn_type == 'PV':
            cnt_pv+=1
        elif mpn_type == 'MF':
            cnt_mf+=1
        elif mpn_type == 'PrePMF':
            cnt_prmf+=1
            
            
    writer = pd.ExcelWriter('{}/{}.xlsx'.format(savepath,cell_folder_list[cc]), engine="xlsxwriter")
    
    df_cell = pd.DataFrame(x_data_cell)
    df_cell.columns = header
    df_cell.set_index('Cell Type', inplace=True)
    df_cell.to_excel(writer, sheet_name="All")
    
    df_cell_n = pd.DataFrame(x_data_cell_n)
    df_cell_n.columns = header_n
    df_cell_n.set_index('Cell Type', inplace=True)
    df_cell_n.to_excel(writer, sheet_name="Normal")
    
    df_cell_et = pd.DataFrame(x_data_cell_et)
    df_cell_et.columns = header_et
    df_cell_et.set_index('Cell Type', inplace=True)
    df_cell_et.to_excel(writer, sheet_name="ET")
    
    df_cell_pv = pd.DataFrame(x_data_cell_pv)
    df_cell_pv.columns = header_pv
    df_cell_pv.set_index('Cell Type', inplace=True)
    df_cell_pv.to_excel(writer, sheet_name="PV")
    
    df_cell_mf = pd.DataFrame(x_data_cell_mf)
    df_cell_mf.columns = header_mf
    df_cell_mf.set_index('Cell Type', inplace=True)
    df_cell_mf.to_excel(writer, sheet_name="MF")
    
    df_cell_prmf = pd.DataFrame(x_data_cell_prmf)
    df_cell_prmf.columns = header_prmf
    df_cell_prmf.set_index('Cell Type', inplace=True)
    df_cell_prmf.to_excel(writer, sheet_name="PrePMF")
    
    writer.close()
