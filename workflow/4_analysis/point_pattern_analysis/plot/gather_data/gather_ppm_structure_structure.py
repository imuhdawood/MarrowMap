import os
import numpy as np
import pandas as pd

xl_path = 'path/to/ppm_result/struct_struct'
key_xl = pd.read_csv('path/to/xenium_keys.csv')
save_path = 'path/to/struct_struct'

if not os.path.exists(save_path):
    os.makedirs(save_path)

structure = ['Bone','Fat','Arteriole','Sinusoid']
structure = np.array(structure).reshape(len(structure),1)
MPN_grp = ['Normal','ET','PV','MF']

for i in range(len(MPN_grp)):
    
    mpn_idx = np.where(np.array(key_xl['Diagnosis_2']) == MPN_grp[i])[0]
    
    xl_list = []
    for ii in range(len(mpn_idx)):
        xl_list.append('{}_R{}'.format(np.array(key_xl['Slide_No'])[mpn_idx[ii]],
                                       np.array(key_xl['Xenium_region'])[mpn_idx[ii]].split('Region')[-1]))
    
    xl_grp = []
    df_cnt = 0
    for ii in range(len(xl_list)):
        try:
            xl_temp = pd.read_excel('{}/{}.xlsx'.format(xl_path,xl_list[ii]))
        except:
            continue
        
        xl_grp.append(xl_list[ii])
        xl_bone = pd.read_excel('{}/{}.xlsx'.format(xl_path,xl_list[ii]),'To Bone')
        xl_fat = pd.read_excel('{}/{}.xlsx'.format(xl_path,xl_list[ii]),'To Fat')
        xl_art = pd.read_excel('{}/{}.xlsx'.format(xl_path,xl_list[ii]),'To Arteriole')
        xl_sinu = pd.read_excel('{}/{}.xlsx'.format(xl_path,xl_list[ii]),'To Sinusoid')
        
        if df_cnt == 0:
            df_temp_b = np.array(xl_bone['beta_1'].rank(method='max')).astype(int)
            df_temp_b = (df_temp_b-min(df_temp_b))/(max(df_temp_b)-min(df_temp_b))
            df_temp_f = np.array(xl_fat['beta_1'].rank(method='max')).astype(int)
            df_temp_f = (df_temp_f-min(df_temp_f))/(max(df_temp_f)-min(df_temp_f))
            df_temp_a = np.array(xl_art['beta_1'].rank(method='max')).astype(int)
            df_temp_a = (df_temp_a-min(df_temp_a))/(max(df_temp_a)-min(df_temp_a))
            df_temp_s = np.array(xl_sinu['beta_1'].rank(method='max')).astype(int)
            df_temp_s = (df_temp_s-min(df_temp_s))/(max(df_temp_s)-min(df_temp_s))
            df_cnt = 1
        else:
            df_ttemp_b = np.array(xl_bone['beta_1'].rank(method='max')).astype(int)
            df_ttemp_b = (df_ttemp_b-min(df_ttemp_b))/(max(df_ttemp_b)-min(df_ttemp_b))
            df_temp_b = np.column_stack((df_temp_b,df_ttemp_b))
            df_ttemp_f = np.array(xl_fat['beta_1'].rank(method='max')).astype(int)
            df_ttemp_f = (df_ttemp_f-min(df_ttemp_f))/(max(df_ttemp_f)-min(df_ttemp_f))
            df_temp_f = np.column_stack((df_temp_f,df_ttemp_f))
            df_ttemp_a = np.array(xl_art['beta_1'].rank(method='max')).astype(int)
            df_ttemp_a = (df_ttemp_a-min(df_ttemp_a))/(max(df_ttemp_a)-min(df_ttemp_a))
            df_temp_a = np.column_stack((df_temp_a,df_ttemp_a))
            df_ttemp_s = np.array(xl_sinu['beta_1'].rank(method='max')).astype(int)
            df_ttemp_s = (df_ttemp_s-min(df_ttemp_s))/(max(df_ttemp_s)-min(df_ttemp_s))
            df_temp_s = np.column_stack((df_temp_s,df_ttemp_s))
        
    df_temp_b = np.concatenate((structure,df_temp_b),axis=1)
    df_b = pd.DataFrame(df_temp_b)
    df_b.columns = np.concatenate((['structure'],xl_grp))
    medians_bone = np.median(df_temp_b[:, 1:].astype(float), axis=1)
    df_b['median'] = medians_bone
    
    df_temp_f = np.concatenate((structure,df_temp_f),axis=1)
    df_f = pd.DataFrame(df_temp_f)
    df_f.columns = np.concatenate((['structure'],xl_grp))
    medians_fat = np.median(df_temp_f[:, 1:].astype(float), axis=1)
    df_f['median'] = medians_fat
    
    df_temp_a = np.concatenate((structure,df_temp_a),axis=1)
    df_a = pd.DataFrame(df_temp_a)
    df_a.columns = np.concatenate((['structure'],xl_grp))
    medians_art = np.median(df_temp_a[:, 1:].astype(float), axis=1)
    df_a['median'] = medians_art
    
    df_temp_s = np.concatenate((structure,df_temp_s),axis=1)
    df_s = pd.DataFrame(df_temp_s)
    df_s.columns = np.concatenate((['structure'],xl_grp))
    medians_sinu = np.median(df_temp_s[:, 1:].astype(float), axis=1)
    df_s['median'] = medians_sinu
    
    writer_strct = pd.ExcelWriter('{}/{}.xlsx'.format(save_path,MPN_grp[i]), engine='xlsxwriter')
    df_b.to_excel(writer_strct, sheet_name="Bone", index=False)
    df_f.to_excel(writer_strct, sheet_name='Fat', index=False)
    df_a.to_excel(writer_strct, sheet_name='Arteriole', index=False)
    df_s.to_excel(writer_strct, sheet_name='Sinusoid', index=False)
    
    writer_strct.close()
