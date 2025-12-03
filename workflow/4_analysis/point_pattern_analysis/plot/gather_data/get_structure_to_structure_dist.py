import os
import pandas as pd
import numpy as np

structure = ['Bone','Fat','Arteriole','Sinusoid']

xl_path = 'metadata/struct_struct'

savepath = 'distance_strct_strct'
if not os.path.exists(savepath):
    os.makedirs(savepath)

xl_list = []
xl_name_list = []

for file in os.listdir(xl_path):
    if file.endswith('csv'):
        xl_list.append(file)
        xl_name_list.append(file.split('.')[0])


for i in range(len(xl_list)):
    xl = pd.read_csv('{}/{}'.format(xl_path,xl_list[i]))
    if len(xl.columns[xl.isna().all()]) > 0:
        continue
    
    writer = pd.ExcelWriter('{}/{}.xlsx'.format(savepath,xl_list[i].split('.csv')[0]), engine="xlsxwriter")
    
    for ii in range(len(structure)):
        col_idx = []
        for iii in range(len(xl.columns)):
            if xl.columns[iii].split('to ')[-1] == structure[ii]:
                col_idx.append(iii)
        
        cnt = 0
        for iii in col_idx:
            from_strct = xl.columns[iii].split(' to')[0]
            dist = np.array(xl.iloc[:,iii])
            dist = dist[~np.isnan(dist)]
            dist = dist.reshape((len(dist),1))
            from_strct_array = np.full(len(dist),from_strct)
            from_strct_array = from_strct_array.reshape((len(from_strct_array),1))
            if cnt == 0:
                to_array = np.concatenate((from_strct_array,dist),axis=1)
                cnt = 1
            else:
                to_array_temp = np.concatenate((from_strct_array,dist),axis=1)
                to_array = np.concatenate((to_array,to_array_temp),axis=0)
        
        df = pd.DataFrame(to_array)
        df.columns = ['From Structure','Distance']
        df.to_excel(writer, sheet_name='{}'.format(structure[ii]),index=False)

    writer.close()
    print('{} out of {} has been done'.format(i+1,len(xl_list)))
