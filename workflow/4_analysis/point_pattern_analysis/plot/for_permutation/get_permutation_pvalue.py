import os
import pandas as pd
import numpy as np
import statistics
import random
from scipy.stats import combine_pvalues

main_data = 'path/to/output'

main_savepth = 'distance_permutation'

compare_type = ['struct_struct', 'struct_celltype', 'struct_cellneighbor',
                'celltype_celltype','cellneighbor_cellneighbor']

for comptype in compare_type:
    if comptype == 'struct_struct':
        structure = ['Bone','Fat','Arteriole','Sinusoid']
        
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
        xl_list = []
        xl_name_list = []

        for file in os.listdir(xl_path):
            if file.endswith('xlsx'):
                xl_list.append(file)
                xl_name_list.append(file.split('.')[0])


        with open('Correspond_MPN.txt') as p:
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
                
                subsavepath = '{}/{}/{}'.format(save_folder,structure[s],MPN_subs[i])
                
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
                    
                    try:
                        x_data = pd.read_excel('{}/{}.xlsx'.format(xl_path,ID_list[mpn_idx[ii]]),structure[s])
                    except:
                        x_data = pd.read_csv('{}/{}.csv'.format(xl_path,ID_list[mpn_idx[ii]]),structure[s])
                    
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
                    
                    print('{}: {}: {} out of {} permutation has been done'.format(comptype, MPN_subs[i], ii+1,len(mpn_idx)))
                    
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
                df.to_csv('{}/pvalue.xlsx'.format(subsavepath),index=False)
                
    elif comptype == 'struct_celltype':
        structure = ['Bone','Fat','Arteriole','Sinusoid']
        
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
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
        
        for s in range(len(structure)):
            for i in range(len(MPN_subs)):
                
                subsavepath = '{}/{}/{}'.format(save_folder,structure[s],MPN_subs[i])
                
                if not os.path.exists(subsavepath):
                    os.makedirs(subsavepath)
                
                mpn_idx = np.where(np.array(MPN_type_list) == MPN_subs[i])[0]
                
            
                for ii in range(len(mpn_idx)):
                    
                    try:
                        x_data = pd.read_excel('{}/{}.xlsx'.format(xl_path,ID_list[mpn_idx[ii]]))
                    except:
                        x_data = pd.read_csv('{}/{}.csv'.format(xl_path,ID_list[mpn_idx[ii]]))
                    
                    x_data_orig = np.zeros((len(cell_type_list),2))
                    x_data_orig = x_data_orig.astype(object)
                    x_data_orig[:,0] = cell_type_list
                    
                    x_data_permut = np.zeros((len(cell_type_list),101))
                    x_data_permut = x_data_permut.astype(object)
                    x_data_permut[:,0] = cell_type_list
                    
                    for iii in range(len(cell_type_list)):
                        type_idx = np.where(np.array(x_data['Cell Type']) == cell_type_list[iii])[0]
                        dist_grp = np.array(x_data['Distance to {}'.format(structure[s])])[type_idx]
                        orig_med_dist = statistics.median(dist_grp)
                        x_data_orig[iii,1] = orig_med_dist
                    
                    for iii in range(0,100):
                        type_temp = np.array(x_data['Cell Type'])
                        random.shuffle(type_temp)
                        for iiii in range(len(cell_type_list)):
                            type_idx = np.where(type_temp == cell_type_list[iiii])[0]
                            bone_dist_grp = np.array(x_data['Distance to Bone'])[type_idx]
                            med_dist_bone = statistics.median(bone_dist_grp)
                            x_data_permut[iiii,iii+1] = med_dist_bone
                    
                    for iii in range(len(x_data_orig)):
                        p_ = len(np.where(np.array(x_data_permut)[iii][1:] < x_data_orig[iii,1])[0])/100
                        pvalue[iii,ii+1] = p_
                    
                    print('{}: {} out of {} permutation has been done'.format(MPN_subs[i], ii+1,len(mpn_idx)))
                
                p_stouff = np.zeros((len(cell_type_list),2))
                p_stouff = p_stouff.astype(object)
                p_stouff[:,0] = cell_type_list
                
                for ii in range(len(pvalue)):
                    stouffer_p = combine_pvalues(np.clip(pvalue[ii][1:].astype(float),1e-16,1-1e-16), 
                                                      method = 'stouffer')[1]
                    
                    p_stouff[ii,1] = stouffer_p
                    
                df = pd.DataFrame(p_stouff)
                df.to_csv('{}/pvalue.csv'.format(subsavepath))
                
    elif comptype == 'struct_cellneighbor':
        
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        
        structure = ['Bone','Fat','Arteriole','Sinusoid']
        cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']
        
        with open('Correspond_MPN.txt') as p:
            strings_p = p.readlines()
            ID_list_orig = [string.split()[0] for string in strings_p]
            MPN_type_list_orig = [string.split()[1] for string in strings_p]
            p.close()
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
        xl_list = []

        for file in os.listdir(xl_path):
            xl_list.append(file)
        
        for s in range(len(structure)):
            ID_list = []
            MPN_type_list = []
            
            for xx in range(len(ID_list_orig)):
                
                try:
                    if '{}.xlsx'.format(ID_list_orig[xx]) in xl_list:
                        ID_list.append(ID_list_orig[xx])
                        MPN_type_list.append(MPN_type_list_orig[xx])
                except:
                    if '{}.csv'.format(ID_list_orig[xx]) in xl_list:
                        ID_list.append(ID_list_orig[xx])
                        MPN_type_list.append(MPN_type_list_orig[xx])
                        
            MPN_subs = np.unique(MPN_type_list)
            
            for i in range(len(MPN_subs)):
                subsavepath = '{}/{}/{}'.format(save_folder,structure[s],MPN_subs[i])
                
                if not os.path.exists(subsavepath):
                    os.makedirs(subsavepath)
                
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
                    
                    try:
                        x_data = pd.read_excel('{}/{}/{}.xlsx'.format(xl_path,structure[s],ID_list[mpn_idx[ii]]))
                    except:
                        x_data = pd.read_csv('{}/{}/{}.csv'.format(xl_path,structure[s],ID_list[mpn_idx[ii]]))
                    
                    dataname = ID_list[mpn_idx[ii]].split('.xlsx')[0]
                    
                    for iii in range(len(cluster_type_list)):
                        type_idx = np.where(np.array(x_data['Cluster Type']) == int(cluster_type_list[iii]))[0]
                        dist_grp = np.array(x_data['Distance to {}'.format(structure[s].capitalize())])[type_idx]
                        # pdb.set_trace()
                        orig_med_dist = statistics.median(dist_grp)
                        x_data_distnc[iii,ii+1] = orig_med_dist
                        x_data_orig[iii,1] = orig_med_dist
                    
                    for iii in range(0,100):
                        type_temp = np.array(x_data['Cluster Type'])
                        random.shuffle(type_temp)
                        for iiii in range(len(cluster_type_list)):
                            type_idx = np.where(type_temp == int(cluster_type_list[iiii]))[0]
                            dist_grp = np.array(x_data['Distance to {}'.format(structure[s].capitalize())])[type_idx]
                            med_dist = statistics.median(dist_grp)
                            x_data_permut[iiii,iii+1] = med_dist
                    
                    for iii in range(len(x_data_orig)):
                        p_ = len(np.where(np.array(x_data_permut)[iii][1:] < x_data_orig[iii,1])[0])/100
                        pvalue[iii,ii+1] = p_
                    
                    print('{} - {}: {} out of {} permutation has been done'.format(structure[s],MPN_subs[i], ii+1,len(mpn_idx)))
                    
                p_stouff = np.zeros((len(cluster_type_list),2))
                p_stouff = p_stouff.astype(object)
                p_stouff[:,0] = cluster_type_list
            
                for ii in range(len(pvalue)):
                    stouffer_p = combine_pvalues(np.clip(pvalue[ii][1:].astype(float),1e-16,1-1e-16), 
                                                      method = 'stouffer')[1]
                    
                    p_stouff[ii,1] = stouffer_p
                    
            
                df = pd.DataFrame(p_stouff)
                df.to_csv('{}/pvalue.csv'.format(subsavepath))
        
    elif comptype == 'celltype_celltype':
        
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
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
                
                subsavepath = '{}/{}/{}'.format(save_folder,cell_type_list[cc],MPN_subs[i])
                
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
                    
                    try:
                        x_data = pd.read_excel('{}/{}/{}.xlsx'.format(xl_path,cell_type_list[cc],ID_list[mpn_idx[ii]]))
                    except:
                        x_data = pd.read_csv('{}/{}/{}.csv'.format(xl_path,cell_type_list[cc],ID_list[mpn_idx[ii]]))
                    
                    dataname = ID_list[mpn_idx[ii]].split('.xlsx')[0]
                    
                    for iii in range(len(cell_type_list)):
                        type_idx = np.where(np.array(x_data['Cell Type']) == cell_type_list[iii])[0]
                        dist_grp = np.array(x_data['Distance to Cell'])[type_idx]
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
                df.to_csv('{}/pvalue.csv'.format(subsavepath))
                
    elif comptype == 'cellneighbor_cellneighbor':
        
        save_folder = '{}/{}'.format(main_savepth,comptype)
        
        with open('Correspond_MPN.txt') as p:
            strings_p = p.readlines()
            ID_list_orig = [string.split()[0] for string in strings_p]
            MPN_type_list_orig = [string.split()[1] for string in strings_p]
            p.close()
        
        xl_path = '{}/{}'.format(main_data,comptype)
        
        xl_list = []

        for file in os.listdir(xl_path):
            xl_list.append(file)
        
        cluster_type_list = ['0','1','2','3','4','5','6','7','8','9']
        
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
                subsavepath = '{}/{}/{}'.format(save_folder,cluster_name,MPN_subs[i])
                
                if not os.path.exists(subsavepath):
                    os.makedirs(subsavepath)
                
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
                    
                    try:
                        x_data = pd.read_excel('{}/{}/{}.xlsx'.format(xl_path,cluster_name,ID_list[mpn_idx[ii]]))
                    except:
                        x_data = pd.read_csv('{}/{}/{}.csv'.format(xl_path,cluster_name,ID_list[mpn_idx[ii]]))
                    
                    dataname = ID_list[mpn_idx[ii]].split('.xlsx')[0]
                    
                    for iii in range(len(cluster_type_list)):
                        type_idx = np.where(np.array(x_data['Cluster Type']) == int(cluster_type_list[iii]))[0]
                        dist_grp = np.array(x_data['Distance to Cluster'])[type_idx]
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
                df.to_csv('{}/pvalue.xlsx'.format(subsavepath))
