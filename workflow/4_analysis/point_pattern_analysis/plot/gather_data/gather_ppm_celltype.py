import os
import pandas as pd
import numpy as np

xl_path = 'path/to/ppm_result/struct_celltype'

savepath = 'path/to/struct_celltype'

if not os.path.exists(savepath):
    os.makedirs(savepath)


with open('xenium_pipeline/celltypes.txt') as p:
    strings_p = p.readlines()
    cell_type_list = [string.split()[0] for string in strings_p]
    p.close()

cell_type_list.remove('CD69')

correspond_mpn = pd.read_excel('xenium_pipeline/Correspond_MPN.xlsx')

mpn_list = np.unique(correspond_mpn['MPN'])

mpn_list = mpn_list[mpn_list != 'PrePMF'] 

struct_grp = ['Bone','Fat','Arteriole','Sinusoid']
# struct_grp = ['Arteriole','Sinusoid']

for g in range(0,len(struct_grp)):
    struct_name = struct_grp[g]
    
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
        x_data_cell = np.zeros((len(cell_type_list),numb_mpn+1))
        x_data_cell = x_data_cell.astype(object)
        x_data_cell[:,0] = cell_type_list
        
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
                type_idx = np.where(np.array(cell_type_list) == xl_whole['Cell Type'][sort_[ii]])[0][0]
                rank_ = ii/max(sort_)
                x_data_cell[type_idx,cnt+1] = rank_
            cnt+=1
            
        df_cell = pd.DataFrame(x_data_cell)
        df_cell = df_cell.rename(columns={df_cell.columns[0]: 'Cell Type'})
        df_cell.set_index('Cell Type', inplace=True)
        df_cell.columns = header
        df_cell.to_excel(writer_, sheet_name='{}'.format(mpn_list[mp]))
            
    writer_.close()
