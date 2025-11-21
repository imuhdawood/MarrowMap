import os
import pandas as pd
import numpy as np
import statistics
import random
from scipy.stats import combine_pvalues

structure = ['Bone','Fat','Arteriole','Sinusoid']

xl_path = 'xenium_pipeline/distance_strct_strct_patchbased_v2'

savepath = 'xenium_pipeline/permut_strct_strct_patchbased_v2_filtered'

if not os.path.exists(savepath):
    os.makedirs(savepath)

xl_list = []
xl_name_list = []

for file in os.listdir(xl_path):
    if file.endswith('xlsx'):
        xl_list.append(file)
        xl_name_list.append(file.split('.')[0])


with open('C:/work/work/code/xenium_pipeline/Correspond_MPN.txt') as p:
    strings_p = p.readlines()
    ID_list_temp = [string.split()[0] for string in strings_p]
    MPN_type_list_temp = [string.split()[1] for string in strings_p]
    p.close()
    

ID_list = []
MPN_type_list = []
for i in range(len(ID_list_temp)):
    if ID_list_temp[i] in xl_name_list:
        ID_list.append(ID_list_temp[i])
        MPN_type_list.append(MPN_type_list_temp[i])

MPN_subs = np.unique(MPN_type_list)

#################################################################
# For each MPN subtypes
for s in range(len(structure)):
    from_strct_grp = structure.copy()
        
    for i in range(len(MPN_subs)):
        
        if MPN_subs[i] == 'PrePMF':
            continue
        
        subsavepath = '{}/{}/{}'.format(savepath,structure[s],MPN_subs[i])
        
        if not os.path.exists(subsavepath):
            os.makedirs(subsavepath)
        
        mpn_idx = np.where(np.array(MPN_type_list) == MPN_subs[i])[0]
        
        pvalue = np.zeros((4,len(mpn_idx)+1))
        pvalue = pvalue.astype(object)
        pvalue[:,0] = from_strct_grp
    
        x_data_distnc = np.zeros((4,len(mpn_idx)+1))
        x_data_distnc = x_data_distnc.astype(object)
        x_data_distnc[:,0] = from_strct_grp
    
        
        for ii in range(len(mpn_idx)):
            x_data_orig = np.zeros((4,2))
            x_data_orig = x_data_orig.astype(object)
            x_data_orig[:,0] = from_strct_grp
            
            x_data_permut = np.zeros((4,101))
            x_data_permut = x_data_permut.astype(object)
            x_data_permut[:,0] = from_strct_grp
            
            x_data = pd.read_excel('{}/{}.xlsx'.format(xl_path,ID_list[mpn_idx[ii]]),structure[s])
            
            dataname = ID_list[mpn_idx[ii]].split('.xlsx')[0]
            
            for iii in range(len(from_strct_grp)):
                type_idx = np.where(np.array(x_data['From Structure']) == from_strct_grp[iii])[0]
                dist_grp = np.array(x_data['Distance'])[type_idx]
                orig_med_dist = statistics.median(dist_grp)
                x_data_distnc[iii,ii+1] = orig_med_dist
                x_data_orig[iii,1] = orig_med_dist
            
            for iii in range(0,100):
                type_temp = np.array(x_data['From Structure'])
                random.shuffle(type_temp)
                for iiii in range(len(from_strct_grp)):
                    type_idx = np.where(type_temp == from_strct_grp[iiii])[0]
                    dist_grp = np.array(x_data['Distance'])[type_idx]
                    med_dist = statistics.median(dist_grp)
                    x_data_permut[iiii,iii+1] = med_dist
                
            for iii in range(len(x_data_orig)):
                p_ = len(np.where(np.array(x_data_permut)[iii][1:] < x_data_orig[iii,1])[0])/100
                pvalue[iii,ii+1] = p_
            
            print('{}: {} out of {} permutation has been done'.format(MPN_subs[i], ii+1,len(mpn_idx)))
            
        p_stouff = np.zeros((4,2))
        p_stouff = p_stouff.astype(object)
        p_stouff[:,0] = from_strct_grp
    
        for ii in range(len(pvalue)):
            stouffer_p = combine_pvalues(np.clip(pvalue[ii][1:].astype(float),1e-16,1-1e-16), 
                                              method = 'stouffer')[1]
            
            
            p_stouff[ii,1] = stouffer_p
    
        df = pd.DataFrame(p_stouff)
        df.columns = ['Structure From','P Value']
        df['Structure To'] = np.full((len(df),1),structure[s])
        df.to_excel('{}/permut.xlsx'.format(subsavepath),index=False)
