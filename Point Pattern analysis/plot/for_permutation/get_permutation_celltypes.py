import os
import pandas as pd
import numpy as np
import statistics
import random
from scipy.stats import combine_pvalues

xl_path = 'path/to/csv' #From get_metadata

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

with open('Correspond_MPN.txt') as p:
    strings_p = p.readlines()
    ID_list = [string.split()[0] for string in strings_p]
    MPN_type_list = [string.split()[1] for string in strings_p]
    p.close()
    
MPN_subs = np.unique(MPN_type_list)

for cc in range(len(cell_type_list)):
    for i in range(len(MPN_subs)):
        
        subsavepath = '{}/{}/{}'.format(savepath,cell_type_list[cc],MPN_subs[i])
        
        if os.path.exists('{}/distance.xlsx'.format(subsavepath)):
            continue
        
        if not os.path.exists(subsavepath):
            os.makedirs(subsavepath)
        
        mpn_idx = np.where(np.array(MPN_type_list) == MPN_subs[i])[0]
        
        pvalue = np.zeros((len(cell_type_list),len(mpn_idx)+1))
        pvalue = pvalue.astype(object)
        pvalue[:,0] = cell_type_list
    
        x_data_distnc = np.zeros((len(cell_type_list),len(mpn_idx)+1))
        x_data_distnc = x_data_distnc.astype(object)
        x_data_distnc[:,0] = cell_type_list
    
        for ii in range(len(mpn_idx)):
            x_data_orig = np.zeros((len(cell_type_list),2))
            x_data_orig = x_data_orig.astype(object)
            x_data_orig[:,0] = cell_type_list
            
            x_data_permut = np.zeros((len(cell_type_list),101))
            x_data_permut = x_data_permut.astype(object)
            x_data_permut[:,0] = cell_type_list
            
            x_data = pd.read_excel('{}/{}/{}.xlsx'.format(xl_path,cell_type_list[cc],ID_list[mpn_idx[ii]]))
            
            dataname = ID_list[mpn_idx[ii]].split('.xlsx')[0]
            
            for iii in range(len(cell_type_list)):
                type_idx = np.where(np.array(x_data['Cell Type']) == cell_type_list[iii])[0]
                dist_grp = np.array(x_data['Distance to Cell'])[type_idx]
                # pdb.set_trace()
                orig_med_dist = statistics.median(dist_grp)
                x_data_distnc[iii,ii+1] = orig_med_dist
                x_data_orig[iii,1] = orig_med_dist
                
            for iii in range(0,100):
                type_temp = np.array(x_data['Cell Type'])
                random.shuffle(type_temp)
                for iiii in range(len(cell_type_list)):
                    type_idx = np.where(type_temp == cell_type_list[iiii])[0]
                    dist_grp = np.array(x_data['Distance to Cell'])[type_idx]
                    med_dist = statistics.median(dist_grp)
                    x_data_permut[iiii,iii+1] = med_dist
                
            for iii in range(len(x_data_orig)):
                p_ = len(np.where(np.array(x_data_permut)[iii][1:] < x_data_orig[iii,1])[0])/100
                pvalue[iii,ii+1] = p_
            
            print('{} - {}: {} out of {} permutation has been done'.format(cell_type_list[cc],MPN_subs[i], ii+1,len(mpn_idx)))
            
        p_stouff = np.zeros((len(cell_type_list),2))
        p_stouff = p_stouff.astype(object)
        p_stouff[:,0] = cell_type_list
    
        for ii in range(len(pvalue)):
            stouffer_p = combine_pvalues(np.clip(pvalue[ii][1:].astype(float),1e-16,1-1e-16), 
                                              method = 'stouffer')[1]
            
            p_stouff[ii,1] = stouffer_p
            
        df = pd.DataFrame(p_stouff)
        df.to_excel('{}/permut.xlsx'.format(subsavepath))
        
        df_dist = pd.DataFrame(x_data_orig)
        df_dist.to_excel('{}/distance.xlsx'.format(subsavepath))
    
        print('{} - {} has been done'.format(cell_type_list[cc],MPN_subs[i]))
        

len_xl_list = len(os.listdir('{}/{}'.format(xl_path,cell_type_list[0])))

xl_grp = []

for i in range(0,len_xl_list):
    xl_grp.append(os.listdir('{}/{}'.format(xl_path,cell_type_list[0]))[i])

for cc in range(len(cell_type_list)):
    pvalue = np.zeros((len(cell_type_list),len_xl_list+1))
    pvalue = pvalue.astype(object)
    pvalue[:,0] = cell_type_list

    x_data_distnc = np.zeros((len(cell_type_list),len_xl_list+1))
    x_data_distnc = x_data_distnc.astype(object)
    x_data_distnc[:,0] = cell_type_list
    subsavepath = '{}/{}/All'.format(savepath,cell_type_list[cc])
    if not os.path.exists(subsavepath):
        os.makedirs(subsavepath)
        
    for i in range(len(xl_grp)):
        x_data_orig = np.zeros((len(cell_type_list),2))
        x_data_orig = x_data_orig.astype(object)
        x_data_orig[:,0] = cell_type_list
        
        x_data_permut = np.zeros((len(cell_type_list),101))
        x_data_permut = x_data_permut.astype(object)
        x_data_permut[:,0] = cell_type_list
        
        dataname = xl_grp[i].split('.xlsx')[0]
        x_data = pd.read_excel('{}/{}/{}.xlsx'.format(xl_path,cell_type_list[cc],dataname))
        for ii in range(len(cell_type_list)):
            type_idx = np.where(np.array(x_data['Cell Type']) == cell_type_list[ii])[0]
            dist_grp = np.array(x_data['Distance to Cell'])[type_idx]
            # pdb.set_trace()
            orig_med_dist = statistics.median(dist_grp)
            x_data_distnc[ii,i+1] = orig_med_dist
            x_data_orig[ii,1] = orig_med_dist
        
        for ii in range(0,100):
            type_temp = np.array(x_data['Cell Type'])
            random.shuffle(type_temp)
            for iii in range(len(cell_type_list)):
                type_idx = np.where(type_temp == cell_type_list[iii])[0]
                dist_grp = np.array(x_data['Distance to Cell'])[type_idx]
                med_dist = statistics.median(dist_grp)
                x_data_permut[iii,ii+1] = med_dist
        
        for ii in range(len(x_data_orig)):
            p_ = len(np.where(np.array(x_data_permut)[ii][1:] < x_data_orig[ii,1])[0])/100
            pvalue[ii,i+1] = p_
        
        print('{} out of {} permutation has been done'.format(i+1,len(xl_list)))

    p_stouff = np.zeros((len(cell_type_list),2))
    p_stouff = p_stouff.astype(object)
    p_stouff[:,0] = cell_type_list
    
    for i in range(len(pvalue)):
        stouffer_p = combine_pvalues(np.clip(pvalue[i][1:].astype(float),1e-16,1-1e-16), 
                                          method = 'stouffer')[1]
        
        p_stouff[i,1] = stouffer_p
        
    # combined_p_value, combined_z_score = stouffers_method(bone_pvalue[0])
    df = pd.DataFrame(p_stouff)
    df.to_excel('{}/permut.xlsx'.format(subsavepath))
    df_dist = pd.DataFrame(x_data_orig)
    df_dist.to_excel('{}/distance.xlsx'.format(subsavepath))
