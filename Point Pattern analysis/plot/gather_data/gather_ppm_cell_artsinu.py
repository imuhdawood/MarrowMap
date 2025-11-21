import os
import pandas as pd
import numpy as np

art_xl_path = 'C:/work/work/code/xenium_pipeline/ppm_results_anno3_stromal_art'
sinu_xl_path = 'C:/work/work/code/xenium_pipeline/ppm_results_anno3_stromal_sinu'

art_savepath = 'C:/work/work/code/xenium_pipeline/gathered_ppm_downsampled_art_anno3_stromal'
sinu_savepath = 'C:/work/work/code/xenium_pipeline/gathered_ppm_downsampled_sinu_anno3_stromal'

if not os.path.exists(art_savepath):
    os.makedirs(art_savepath)

if not os.path.exists(sinu_savepath):
    os.makedirs(sinu_savepath)

art_xl_list = []

for file in os.listdir(art_xl_path):
    art_xl_list.append(file)

sinu_xl_list = []

for file in os.listdir(sinu_xl_path):
    sinu_xl_list.append(file)

with open('celltypes.txt') as p:
    strings_p = p.readlines()
    cell_type_list = [string.split()[0] for string in strings_p]
    p.close()

cell_type_list.remove('CD69')

correspond_mpn = pd.read_excel('Correspond_MPN.xlsx')

correspond_mpn_art = np.zeros((len(art_xl_list),2))
correspond_mpn_art = correspond_mpn_art.astype(object)
for i in range(len(correspond_mpn_art)):
    correspond_mpn_art[i,0] = art_xl_list[i].split('.')[0]
    idx = np.where(np.array(correspond_mpn) == art_xl_list[i].split('.')[0])[0][0]
    correspond_mpn_art[i,1] = correspond_mpn['MPN'][idx]

header_art = ['ID','MPN']

correspond_mpn_art = pd.DataFrame(correspond_mpn_art)
correspond_mpn_art.columns = header_art

mpn_list = np.unique(correspond_mpn_art['MPN'])

##############################################################################
#Arteriole

x_data_cell_art = np.zeros((len(cell_type_list),len(art_xl_list)+1))
x_data_cell_art = x_data_cell_art.astype(object)
x_data_cell_art[:,0] = cell_type_list

for mp in mpn_list:
    numb_mpn = len(np.where(np.array(correspond_mpn_art['MPN']) == mp)[0])
    if mp =='Normal':
        x_data_cell_n_art = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_n_art = x_data_cell_n_art.astype(object)
        x_data_cell_n_art[:,0] = cell_type_list
    elif mp =='ET':
        x_data_cell_et_art = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_et_art = x_data_cell_et_art.astype(object)
        x_data_cell_et_art[:,0] = cell_type_list
    elif mp =='PV':
        x_data_cell_pv_art = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_pv_art = x_data_cell_pv_art.astype(object)
        x_data_cell_pv_art[:,0] = cell_type_list
    elif mp =='MF':
        x_data_cell_mf_art = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_mf_art = x_data_cell_mf_art.astype(object)
        x_data_cell_mf_art[:,0] = cell_type_list
    # elif mp =='PrePMF':
    #     x_data_cell_prmf_art = np.zeros((len(cell_type_list),numb_mpn+1))
    #     x_data_cell_prmf_art = x_data_cell_prmf_art.astype(object)
    #     x_data_cell_prmf_art[:,0] = cell_type_list

header = ['Cell Type']
header_n = ['Cell Type']
header_et = ['Cell Type']
header_pv = ['Cell Type']
header_mf = ['Cell Type']
# header_prmf = ['Cell Type']

cnt_n_art = 0
cnt_et_art = 0
cnt_pv_art = 0
cnt_mf_art = 0
cnt_prmf_art = 0

cnt_n_sinu = 0
cnt_et_sinu = 0
cnt_pv_sinu = 0
cnt_mf_sinu = 0
cnt_prmf_sinu = 0

for i in range(len(art_xl_list)):
    dataname = art_xl_list[i].split('.csv')[0]
    mpn_type = correspond_mpn_art['MPN'][np.where(np.array(correspond_mpn_art['ID'] == dataname))[0][0]]
    header.append('{} {}'.format(dataname,mpn_type))
    
    if mpn_type == 'Normal':
        header_n.append('{} {}'.format(dataname,mpn_type))
    elif mpn_type == 'ET':
        header_et.append('{} {}'.format(dataname,mpn_type))
    elif mpn_type == 'PV':
        header_pv.append('{} {}'.format(dataname,mpn_type))
    elif mpn_type == 'MF':
        header_mf.append('{} {}'.format(dataname,mpn_type))
    
    xl_whole = pd.read_excel('{}/{}'.format(art_xl_path,art_xl_list[i]))
    
    if len(np.where(np.array(xl_whole['name'] == 'CD69'))[0]) > 0:
        c_idx = np.where(np.array(xl_whole['name'] == 'CD69'))[0][0]
        xl_whole = xl_whole.drop([c_idx])
        xl_whole = xl_whole.reset_index(drop=True)
    
    art_sort = np.argsort(xl_whole['distance_to_arteriole'])
    
    for ii in range(len(art_sort)):
        type_idx = np.where(np.array(cell_type_list) == xl_whole['name'][art_sort[ii]])[0][0]
        art_rank = ii/max(art_sort)
        x_data_cell_art[type_idx,i+1] = art_rank
        
        if mpn_type == 'Normal':
            x_data_cell_n_art[type_idx,cnt_n_art+1] = art_rank
        elif mpn_type == 'ET':
            x_data_cell_et_art[type_idx,cnt_et_art+1] = art_rank
        elif mpn_type == 'PV':
            x_data_cell_pv_art[type_idx,cnt_pv_art+1] = art_rank
        elif mpn_type == 'MF':
            x_data_cell_mf_art[type_idx,cnt_mf_art+1] = art_rank
        
    
    if mpn_type == 'Normal':
        cnt_n_art+=1
    elif mpn_type == 'ET':
        cnt_et_art+=1
    elif mpn_type == 'PV':
        cnt_pv_art+=1
    elif mpn_type == 'MF':
        cnt_mf_art+=1
        
writer_art = pd.ExcelWriter('{}/arteriole.xlsx'.format(art_savepath), engine="xlsxwriter")

df_cell_art = pd.DataFrame(x_data_cell_art)
df_cell_art.columns = header
df_cell_art.set_index('Cell Type', inplace=True)
df_cell_art.to_excel(writer_art, sheet_name="All")
# df_cell.to_excel('{}/{}.xlsx'.format(savepath,cell_folder_list[cc]))

df_cell_n_art = pd.DataFrame(x_data_cell_n_art)
df_cell_n_art.columns = header_n
df_cell_n_art.set_index('Cell Type', inplace=True)
df_cell_n_art.to_excel(writer_art, sheet_name="Normal")

df_cell_et_art = pd.DataFrame(x_data_cell_et_art)
df_cell_et_art.columns = header_et
df_cell_et_art.set_index('Cell Type', inplace=True)
df_cell_et_art.to_excel(writer_art, sheet_name="ET")

df_cell_pv_art = pd.DataFrame(x_data_cell_pv_art)
df_cell_pv_art.columns = header_pv
df_cell_pv_art.set_index('Cell Type', inplace=True)
df_cell_pv_art.to_excel(writer_art, sheet_name="PV")

df_cell_mf_art = pd.DataFrame(x_data_cell_mf_art)
df_cell_mf_art.columns = header_mf
df_cell_mf_art.set_index('Cell Type', inplace=True)
df_cell_mf_art.to_excel(writer_art, sheet_name="MF")

writer_art.close()

#############################################################################

#Sinusoid

correspond_mpn_sinu = np.zeros((len(sinu_xl_list),2))
correspond_mpn_sinu = correspond_mpn_sinu.astype(object)
for i in range(len(correspond_mpn_sinu)):
    correspond_mpn_sinu[i,0] = sinu_xl_list[i].split('.')[0]
    idx = np.where(np.array(correspond_mpn) == sinu_xl_list[i].split('.')[0])[0][0]
    correspond_mpn_sinu[i,1] = correspond_mpn['MPN'][idx]

header_sinu = ['ID','MPN']

correspond_mpn_sinu = pd.DataFrame(correspond_mpn_sinu)
correspond_mpn_sinu.columns = header_sinu

x_data_cell_sinu = np.zeros((len(cell_type_list),len(sinu_xl_list)+1))
x_data_cell_sinu = x_data_cell_sinu.astype(object)
x_data_cell_sinu[:,0] = cell_type_list

for mp in mpn_list:
    numb_mpn = len(np.where(np.array(correspond_mpn_sinu['MPN']) == mp)[0])
    if mp =='Normal':
        x_data_cell_n_sinu = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_n_sinu = x_data_cell_n_sinu.astype(object)
        x_data_cell_n_sinu[:,0] = cell_type_list
    elif mp =='ET':
        x_data_cell_et_sinu = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_et_sinu = x_data_cell_et_sinu.astype(object)
        x_data_cell_et_sinu[:,0] = cell_type_list
    elif mp =='PV':
        x_data_cell_pv_sinu = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_pv_sinu = x_data_cell_pv_sinu.astype(object)
        x_data_cell_pv_sinu[:,0] = cell_type_list
    elif mp =='MF':
        x_data_cell_mf_sinu = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_mf_sinu = x_data_cell_mf_sinu.astype(object)
        x_data_cell_mf_sinu[:,0] = cell_type_list
    
header = ['Cell Type']
header_n = ['Cell Type']
header_et = ['Cell Type']
header_pv = ['Cell Type']
header_mf = ['Cell Type']
header_prmf = ['Cell Type']

cnt_n_sinu = 0
cnt_et_sinu = 0
cnt_pv_sinu = 0
cnt_mf_sinu = 0
cnt_prmf_sinu = 0

for i in range(len(sinu_xl_list)):
    dataname = sinu_xl_list[i].split('.csv')[0]
    mpn_type = correspond_mpn_sinu['MPN'][np.where(np.array(correspond_mpn_sinu['ID'] == dataname))[0][0]]
    header.append('{} {}'.format(dataname,mpn_type))
    
    if mpn_type == 'Normal':
        header_n.append('{} {}'.format(dataname,mpn_type))
    elif mpn_type == 'ET':
        header_et.append('{} {}'.format(dataname,mpn_type))
    elif mpn_type == 'PV':
        header_pv.append('{} {}'.format(dataname,mpn_type))
    elif mpn_type == 'MF':
        header_mf.append('{} {}'.format(dataname,mpn_type))
    
    xl_whole = pd.read_excel('{}/{}'.format(sinu_xl_path,sinu_xl_list[i]))
    
    if len(np.where(np.array(xl_whole['name'] == 'CD69'))[0]) > 0:
        c_idx = np.where(np.array(xl_whole['name'] == 'CD69'))[0][0]
        xl_whole = xl_whole.drop([c_idx])
        xl_whole = xl_whole.reset_index(drop=True)
    
    sinu_sort = np.argsort(xl_whole['distance_to_sinusoid'])
    
    for ii in range(len(sinu_sort)):
        type_idx = np.where(np.array(cell_type_list) == xl_whole['name'][sinu_sort[ii]])[0][0]
        sinu_rank = ii/max(sinu_sort)
        x_data_cell_sinu[type_idx,i+1] = sinu_rank
        
        if mpn_type == 'Normal':
            x_data_cell_n_sinu[type_idx,cnt_n_sinu+1] = sinu_rank
        elif mpn_type == 'ET':
            x_data_cell_et_sinu[type_idx,cnt_et_sinu+1] = sinu_rank
        elif mpn_type == 'PV':
            x_data_cell_pv_sinu[type_idx,cnt_pv_sinu+1] = sinu_rank
        elif mpn_type == 'MF':
            x_data_cell_mf_sinu[type_idx,cnt_mf_sinu+1] = sinu_rank
        
    if mpn_type == 'Normal':
        cnt_n_sinu+=1
    elif mpn_type == 'ET':
        cnt_et_sinu+=1
    elif mpn_type == 'PV':
        cnt_pv_sinu+=1
    elif mpn_type == 'MF':
        cnt_mf_sinu+=1
    
writer_sinu = pd.ExcelWriter('{}/sinusoid.xlsx'.format(sinu_savepath), engine="xlsxwriter")

df_cell_sinu = pd.DataFrame(x_data_cell_sinu)
df_cell_sinu.columns = header
df_cell_sinu.set_index('Cell Type', inplace=True)
df_cell_sinu.to_excel(writer_sinu, sheet_name="All")

df_cell_n_sinu = pd.DataFrame(x_data_cell_n_sinu)
df_cell_n_sinu.columns = header_n
df_cell_n_sinu.set_index('Cell Type', inplace=True)
df_cell_n_sinu.to_excel(writer_sinu, sheet_name="Normal")

df_cell_et_sinu = pd.DataFrame(x_data_cell_et_sinu)
df_cell_et_sinu.columns = header_et
df_cell_et_sinu.set_index('Cell Type', inplace=True)
df_cell_et_sinu.to_excel(writer_sinu, sheet_name="ET")

df_cell_pv_sinu = pd.DataFrame(x_data_cell_pv_sinu)
df_cell_pv_sinu.columns = header_pv
df_cell_pv_sinu.set_index('Cell Type', inplace=True)
df_cell_pv_sinu.to_excel(writer_sinu, sheet_name="PV")

df_cell_mf_sinu = pd.DataFrame(x_data_cell_mf_sinu)
df_cell_mf_sinu.columns = header_mf
df_cell_mf_sinu.set_index('Cell Type', inplace=True)
df_cell_mf_sinu.to_excel(writer_sinu, sheet_name="MF")

writer_sinu.close()
