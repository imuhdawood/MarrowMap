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

cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']

with open('Correspond_MPN.txt') as p:
    strings_p = p.readlines()
    ID_list_orig = [string.split()[0] for string in strings_p]
    MPN_type_list_orig = [string.split()[1] for string in strings_p]
    p.close()

for s in range(len(cluster_type_list)):
    cluster_name = cluster_type_list[s]
    sub_xlpath = '{}/{}'.format(xl_path,cluster_name)
    xl_list = os.listdir(sub_xlpath)
    
    ID_list = []
    MPN_type_list = []
    
    for xx in range(len(ID_list_orig)):
        if '{}.xlsx'.format(ID_list_orig[xx]) in xl_list:
            ID_list.append(ID_list_orig[xx])
            MPN_type_list.append(MPN_type_list_orig[xx])
            
    MPN_subs = np.unique(MPN_type_list)
    
    for i in range(len(MPN_subs)):
        subsavepath = '{}/{}/{}'.format(savepath,cluster_name,MPN_subs[i])
        
        if not os.path.exists(subsavepath):
            os.makedirs(subsavepath)
        
        if os.path.exists('{}/distance.xlsx'.format(subsavepath)):
            continue
        
        mpn_idx = np.where(np.array(MPN_type_list) == MPN_subs[i])[0]
        
        pvalue = np.zeros((len(cluster_type_list),len(mpn_idx)+1))
        pvalue = pvalue.astype(object)
        pvalue[:,0] = cluster_type_list
    
        x_data_distnc = np.zeros((len(cluster_type_list),len(mpn_idx)+1))
        x_data_distnc = x_data_distnc.astype(object)
        x_data_distnc[:,0] = cluster_type_list
    
        for ii in range(len(mpn_idx)):
            x_data_orig = np.zeros((len(cluster_type_list),2))
            x_data_orig = x_data_orig.astype(object)
            x_data_orig[:,0] = cluster_type_list
            
            x_data_permut = np.zeros((len(cluster_type_list),101))
            x_data_permut = x_data_permut.astype(object)
            x_data_permut[:,0] = cluster_type_list
            
            x_data = pd.read_excel('{}/{}/{}.xlsx'.format(xl_path,cluster_name,ID_list[mpn_idx[ii]]))
            
            dataname = ID_list[mpn_idx[ii]].split('.xlsx')[0]
            
            for iii in range(len(cluster_type_list)):
                type_idx = np.where(np.array(x_data['Cluster Type']) == int(cluster_type_list[iii]))[0]
                dist_grp = np.array(x_data['Distance to Cluster'])[type_idx]
                # pdb.set_trace()
                orig_med_dist = statistics.median(dist_grp)
                x_data_distnc[iii,ii+1] = orig_med_dist
                x_data_orig[iii,1] = orig_med_dist
                
            for iii in range(0,100):
                type_temp = np.array(x_data['Cluster Type'])
                random.shuffle(type_temp)
                for iiii in range(len(cluster_type_list)):
                    type_idx = np.where(type_temp == int(cluster_type_list[iiii]))[0]
                    dist_grp = np.array(x_data['Distance to Cluster'])[type_idx]
                    med_dist = statistics.median(dist_grp)
                    x_data_permut[iiii,iii+1] = med_dist
            
            for iii in range(len(x_data_orig)):
                p_ = len(np.where(np.array(x_data_permut)[iii][1:] < x_data_orig[iii,1])[0])/100
                pvalue[iii,ii+1] = p_
            
            print('{} - {}: {} out of {} permutation has been done'.format(cluster_name,MPN_subs[i], ii+1,len(mpn_idx)))
            
        p_stouff = np.zeros((len(cluster_type_list),2))
        p_stouff = p_stouff.astype(object)
        p_stouff[:,0] = cluster_type_list
    
        for ii in range(len(pvalue)):
            stouffer_p = combine_pvalues(np.clip(pvalue[ii][1:].astype(float),1e-16,1-1e-16), 
                                              method = 'stouffer')[1]
            
            
            p_stouff[ii,1] = stouffer_p
            
        
        df = pd.DataFrame(p_stouff)
        df.to_excel('{}/permut.xlsx'.format(subsavepath))
        
        df_dist = pd.DataFrame(x_data_orig)
        df_dist.to_excel('{}/distance.xlsx'.format(subsavepath))
