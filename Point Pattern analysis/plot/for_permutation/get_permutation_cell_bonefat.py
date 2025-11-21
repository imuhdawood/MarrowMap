import os
import pandas as pd
import numpy as np
import statistics
import random
from scipy.stats import combine_pvalues

xl_path = 'path/to/data' # From get_metadata

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

x_data_bone_orig = np.zeros((len(cell_type_list),6))
x_data_bone_orig = x_data_bone_orig.astype(object)
x_data_bone_orig[:,0] = cell_type_list

#################################################################
# For each MPN subtypes
for i in range(len(MPN_subs)):
    
    subsavepath = '{}/{}'.format(savepath,MPN_subs[i])
    
    if os.path.exists('{}/fat_distance.xlsx'.format(subsavepath)):
        continue
    
    if not os.path.exists(subsavepath):
        os.makedirs(subsavepath)
    
    
    mpn_idx = np.where(np.array(MPN_type_list) == MPN_subs[i])[0]
    
    bone_pvalue = np.zeros((len(cell_type_list),len(mpn_idx)+1))
    bone_pvalue = bone_pvalue.astype(object)
    bone_pvalue[:,0] = cell_type_list

    fat_pvalue = np.zeros((len(cell_type_list),len(mpn_idx)+1))
    fat_pvalue = fat_pvalue.astype(object)
    fat_pvalue[:,0] = cell_type_list
    
    x_data_bone_distnc = np.zeros((len(cell_type_list),len(mpn_idx)+1))
    x_data_bone_distnc = x_data_bone_distnc.astype(object)
    x_data_bone_distnc[:,0] = cell_type_list

    x_data_fat_distnc = np.zeros((len(cell_type_list),len(mpn_idx)+1))
    x_data_fat_distnc = x_data_fat_distnc.astype(object)
    x_data_fat_distnc[:,0] = cell_type_list
    
    for ii in range(len(mpn_idx)):
        x_data_bone_orig = np.zeros((len(cell_type_list),2))
        x_data_bone_orig = x_data_bone_orig.astype(object)
        x_data_bone_orig[:,0] = cell_type_list
        
        x_data_bone_permut = np.zeros((len(cell_type_list),101))
        x_data_bone_permut = x_data_bone_permut.astype(object)
        x_data_bone_permut[:,0] = cell_type_list
        
        x_data_fat_orig = np.zeros((len(cell_type_list),2))
        x_data_fat_orig = x_data_fat_orig.astype(object)
        x_data_fat_orig[:,0] = cell_type_list
        
        x_data_fat_permut = np.zeros((len(cell_type_list),101))
        x_data_fat_permut = x_data_fat_permut.astype(object)
        x_data_fat_permut[:,0] = cell_type_list
        
        x_data = pd.read_excel('{}/{}.xlsx'.format(xl_path,ID_list[mpn_idx[ii]]))
        
        dataname = ID_list[mpn_idx[ii]].split('.xlsx')[0]
        for iii in range(len(cell_type_list)):
            type_idx = np.where(np.array(x_data['Cell Type']) == cell_type_list[iii])[0]
            bone_dist_grp = np.array(x_data['Distance to Bone'])[type_idx]
            # pdb.set_trace()
            orig_med_dist_bone = statistics.median(bone_dist_grp)
            x_data_bone_distnc[iii,ii+1] = orig_med_dist_bone
            x_data_bone_orig[iii,1] = orig_med_dist_bone
            fat_dist_grp = np.array(x_data['Distance to Fat'])[type_idx]
            orig_med_dist_fat = statistics.median(fat_dist_grp)
            x_data_fat_distnc[iii,ii+1] = orig_med_dist_fat
            x_data_fat_orig[iii,1] = orig_med_dist_fat
        
        
        for iii in range(0,100):
            type_temp = np.array(x_data['Cell Type'])
            random.shuffle(type_temp)
            for iiii in range(len(cell_type_list)):
                type_idx = np.where(type_temp == cell_type_list[iiii])[0]
                bone_dist_grp = np.array(x_data['Distance to Bone'])[type_idx]
                med_dist_bone = statistics.median(bone_dist_grp)
                x_data_bone_permut[iiii,iii+1] = med_dist_bone
            
            type_temp = np.array(x_data['Cell Type'])
            random.shuffle(type_temp)
            for iiii in range(len(cell_type_list)):
                type_idx = np.where(type_temp == cell_type_list[iiii])[0]
                fat_dist_grp = np.array(x_data['Distance to Fat'])[type_idx]
                med_dist_fat = statistics.median(fat_dist_grp)
                x_data_fat_permut[iiii,iii+1] = med_dist_fat
        
        
        for iii in range(len(x_data_bone_orig)):
            bone_p = len(np.where(np.array(x_data_bone_permut)[iii][1:] < x_data_bone_orig[iii,1])[0])/100
            bone_pvalue[iii,ii+1] = bone_p
        
        for iii in range(len(x_data_fat_orig)):
            fat_p = len(np.where(np.array(x_data_fat_permut)[iii][1:] < x_data_fat_orig[iii,1])[0])/100
            fat_pvalue[iii,ii+1] = fat_p

        print('{}: {} out of {} permutation has been done'.format(MPN_subs[i], ii+1,len(mpn_idx)))

    p_stouff_bone = np.zeros((len(cell_type_list),2))
    p_stouff_bone = p_stouff_bone.astype(object)
    p_stouff_bone[:,0] = cell_type_list

    p_stouff_fat = np.zeros((len(cell_type_list),2))
    p_stouff_fat = p_stouff_fat.astype(object)
    p_stouff_fat[:,0] = cell_type_list

    for ii in range(len(bone_pvalue)):
        stouffer_p_bone = combine_pvalues(np.clip(bone_pvalue[ii][1:].astype(float),1e-16,1-1e-16), 
                                          method = 'stouffer')[1]
        
        stouffer_p_fat = combine_pvalues(np.clip(fat_pvalue[ii][1:].astype(float),1e-16,1-1e-16), 
                                          method = 'stouffer')[1]
        
        p_stouff_bone[ii,1] = stouffer_p_bone
        p_stouff_fat[ii,1] = stouffer_p_fat
        
    df_bone = pd.DataFrame(p_stouff_bone)
    df_bone.to_excel('{}/bone.xlsx'.format(subsavepath))
    df_fat = pd.DataFrame(p_stouff_fat)
    df_fat.to_excel('{}/fat.xlsx'.format(subsavepath))

    df_dist_bone = pd.DataFrame(x_data_bone_orig)
    df_dist_bone.to_excel('{}/bone_distance.xlsx'.format(subsavepath))

    df_dist_fat = pd.DataFrame(x_data_fat_orig)
    df_dist_fat.to_excel('{}/fat_distance.xlsx'.format(subsavepath))
