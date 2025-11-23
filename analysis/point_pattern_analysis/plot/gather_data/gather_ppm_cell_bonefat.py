import os
import pandas as pd
import numpy as np

xl_path = 'path/to/data'

savepath = 'path/to/output'

if not os.path.exists(savepath):
    os.makedirs(savepath)

xl_list = []

for file in os.listdir(xl_path):
    xl_list.append(file)

with open('celltypes.txt') as p:
    strings_p = p.readlines()
    cell_type_list = [string.split()[0] for string in strings_p]
    p.close()

cell_type_list.remove('CD69')

correspond_mpn = pd.read_excel('Correspond_MPN.xlsx')

mpn_list = np.unique(correspond_mpn['MPN'])

xl_list = []

for file in os.listdir(xl_path):
    xl_list.append(file)

x_data_cell_bone = np.zeros((len(cell_type_list),len(xl_list)+1))
x_data_cell_bone = x_data_cell_bone.astype(object)
x_data_cell_bone[:,0] = cell_type_list

x_data_cell_fat = np.zeros((len(cell_type_list),len(xl_list)+1))
x_data_cell_fat = x_data_cell_fat.astype(object)
x_data_cell_fat[:,0] = cell_type_list

for mp in mpn_list:
    numb_mpn = len(np.where(np.array(correspond_mpn['MPN']) == mp)[0])
    if mp =='Normal':
        x_data_cell_n_bone = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_n_bone = x_data_cell_n_bone.astype(object)
        x_data_cell_n_bone[:,0] = cell_type_list
        x_data_cell_n_fat = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_n_fat = x_data_cell_n_fat.astype(object)
        x_data_cell_n_fat[:,0] = cell_type_list
    elif mp =='ET':
        x_data_cell_et_bone = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_et_bone = x_data_cell_et_bone.astype(object)
        x_data_cell_et_bone[:,0] = cell_type_list
        x_data_cell_et_fat = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_et_fat = x_data_cell_et_fat.astype(object)
        x_data_cell_et_fat[:,0] = cell_type_list
    elif mp =='PV':
        x_data_cell_pv_bone = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_pv_bone = x_data_cell_pv_bone.astype(object)
        x_data_cell_pv_bone[:,0] = cell_type_list
        x_data_cell_pv_fat = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_pv_fat = x_data_cell_pv_fat.astype(object)
        x_data_cell_pv_fat[:,0] = cell_type_list
    elif mp =='MF':
        x_data_cell_mf_bone = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_mf_bone = x_data_cell_mf_bone.astype(object)
        x_data_cell_mf_bone[:,0] = cell_type_list
        x_data_cell_mf_fat = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_mf_fat = x_data_cell_mf_fat.astype(object)
        x_data_cell_mf_fat[:,0] = cell_type_list
    elif mp =='PrePMF':
        x_data_cell_prmf_bone = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_prmf_bone = x_data_cell_prmf_bone.astype(object)
        x_data_cell_prmf_bone[:,0] = cell_type_list
        x_data_cell_prmf_fat = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell_prmf_fat = x_data_cell_prmf_fat.astype(object)
        x_data_cell_prmf_fat[:,0] = cell_type_list

header = ['Cell Type']
header_n = ['Cell Type']
header_et = ['Cell Type']
header_pv = ['Cell Type']
header_mf = ['Cell Type']
header_prmf = ['Cell Type']

cnt_n_bone = 0
cnt_et_bone = 0
cnt_pv_bone = 0
cnt_mf_bone = 0
cnt_prmf_bone = 0

cnt_n_fat = 0
cnt_et_fat = 0
cnt_pv_fat = 0
cnt_mf_fat = 0
cnt_prmf_fat = 0


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
    
    xl_whole = pd.read_excel('{}/{}'.format(xl_path,xl_list[i]))
    
    if len(np.where(np.array(xl_whole['name'] == 'CD69'))[0]) > 0:
        c_idx = np.where(np.array(xl_whole['name'] == 'CD69'))[0][0]
        xl_whole = xl_whole.drop([c_idx])
        xl_whole = xl_whole.reset_index(drop=True)
    
    bone_sort = np.argsort(xl_whole['distance_to_bone'])
    
    fat_sort = np.argsort(xl_whole['distance_to_fat'])
    
    # pdb.set_trace()
    for ii in range(len(bone_sort)):
        type_idx = np.where(np.array(cell_type_list) == xl_whole['name'][bone_sort[ii]])[0][0]
        bone_rank = ii/max(bone_sort)
        x_data_cell_bone[type_idx,i+1] = bone_rank
        
        if mpn_type == 'Normal':
            x_data_cell_n_bone[type_idx,cnt_n_bone+1] = bone_rank
        elif mpn_type == 'ET':
            x_data_cell_et_bone[type_idx,cnt_et_bone+1] = bone_rank
        elif mpn_type == 'PV':
            x_data_cell_pv_bone[type_idx,cnt_pv_bone+1] = bone_rank
        elif mpn_type == 'MF':
            x_data_cell_mf_bone[type_idx,cnt_mf_bone+1] = bone_rank
        elif mpn_type == 'PrePMF':
            x_data_cell_prmf_bone[type_idx,cnt_prmf_bone+1] = bone_rank
    
    
    for ii in range(len(fat_sort)):
        type_idx = np.where(np.array(cell_type_list) == xl_whole['name'][fat_sort[ii]])[0][0]
        fat_rank = ii/max(fat_sort)
        x_data_cell_fat[type_idx,i+1] = fat_rank
        
        if mpn_type == 'Normal':
            x_data_cell_n_fat[type_idx,cnt_n_fat+1] = fat_rank
        elif mpn_type == 'ET':
            x_data_cell_et_fat[type_idx,cnt_et_fat+1] = fat_rank
        elif mpn_type == 'PV':
            x_data_cell_pv_fat[type_idx,cnt_pv_fat+1] = fat_rank
        elif mpn_type == 'MF':
            x_data_cell_mf_fat[type_idx,cnt_mf_fat+1] = fat_rank
        elif mpn_type == 'PrePMF':
            x_data_cell_prmf_fat[type_idx,cnt_prmf_fat+1] = fat_rank
            
    if mpn_type == 'Normal':
        cnt_n_bone+=1
        cnt_n_fat+=1
    elif mpn_type == 'ET':
        cnt_et_bone+=1
        cnt_et_fat+=1
    elif mpn_type == 'PV':
        cnt_pv_bone+=1
        cnt_pv_fat+=1
    elif mpn_type == 'MF':
        cnt_mf_bone+=1
        cnt_mf_fat+=1
    elif mpn_type == 'PrePMF':
        cnt_prmf_bone+=1
        cnt_prmf_fat+=1
        
writer_bone = pd.ExcelWriter('{}/bone.xlsx'.format(savepath), engine="xlsxwriter")

df_cell_bone = pd.DataFrame(x_data_cell_bone)
df_cell_bone.columns = header
df_cell_bone.set_index('Cell Type', inplace=True)
df_cell_bone.to_excel(writer_bone, sheet_name="All")

df_cell_n_bone = pd.DataFrame(x_data_cell_n_bone)
df_cell_n_bone.columns = header_n
df_cell_n_bone.set_index('Cell Type', inplace=True)
df_cell_n_bone.to_excel(writer_bone, sheet_name="Normal")

df_cell_et_bone = pd.DataFrame(x_data_cell_et_bone)
df_cell_et_bone.columns = header_et
df_cell_et_bone.set_index('Cell Type', inplace=True)
df_cell_et_bone.to_excel(writer_bone, sheet_name="ET")

df_cell_pv_bone = pd.DataFrame(x_data_cell_pv_bone)
df_cell_pv_bone.columns = header_pv
df_cell_pv_bone.set_index('Cell Type', inplace=True)
df_cell_pv_bone.to_excel(writer_bone, sheet_name="PV")

df_cell_mf_bone = pd.DataFrame(x_data_cell_mf_bone)
df_cell_mf_bone.columns = header_mf
df_cell_mf_bone.set_index('Cell Type', inplace=True)
df_cell_mf_bone.to_excel(writer_bone, sheet_name="MF")

df_cell_prmf_bone = pd.DataFrame(x_data_cell_prmf_bone)
df_cell_prmf_bone.columns = header_prmf
df_cell_prmf_bone.set_index('Cell Type', inplace=True)
df_cell_prmf_bone.to_excel(writer_bone, sheet_name="PrePMF")

writer_bone.close()

writer_fat = pd.ExcelWriter('{}/fat.xlsx'.format(savepath), engine="xlsxwriter")

df_cell_fat = pd.DataFrame(x_data_cell_fat)
df_cell_fat.columns = header
df_cell_fat.set_index('Cell Type', inplace=True)
df_cell_fat.to_excel(writer_fat, sheet_name="All")

df_cell_n_fat = pd.DataFrame(x_data_cell_n_fat)
df_cell_n_fat.columns = header_n
df_cell_n_fat.set_index('Cell Type', inplace=True)
df_cell_n_fat.to_excel(writer_fat, sheet_name="Normal")

df_cell_et_fat = pd.DataFrame(x_data_cell_et_fat)
df_cell_et_fat.columns = header_et
df_cell_et_fat.set_index('Cell Type', inplace=True)
df_cell_et_fat.to_excel(writer_fat, sheet_name="ET")

df_cell_pv_fat = pd.DataFrame(x_data_cell_pv_fat)
df_cell_pv_fat.columns = header_pv
df_cell_pv_fat.set_index('Cell Type', inplace=True)
df_cell_pv_fat.to_excel(writer_fat, sheet_name="PV")

df_cell_mf_fat = pd.DataFrame(x_data_cell_mf_fat)
df_cell_mf_fat.columns = header_mf
df_cell_mf_fat.set_index('Cell Type', inplace=True)
df_cell_mf_fat.to_excel(writer_fat, sheet_name="MF")

df_cell_prmf_fat = pd.DataFrame(x_data_cell_prmf_fat)
df_cell_prmf_fat.columns = header_prmf
df_cell_prmf_fat.set_index('Cell Type', inplace=True)
df_cell_prmf_fat.to_excel(writer_fat, sheet_name="PrePMF")

writer_fat.close()
