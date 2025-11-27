import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon
from scipy.stats import ranksums
from matplotlib import pyplot as plt
import pdb

type_strct = ['arteriole', 'bone','fat','sinusoid']

xl_path = 'xenium_pipeline/gathered_ppm_anno3_cellneighbor'

# xl_path = 'xenium_pipeline/gathered_ppm_downsampled_cell/Granulocyte_mast.xlsx'

compare_with = ['MF','ET','PV']


for i in range(len(compare_with)):
    save_path = 'xenium_pipeline/delta_norm{}_cellneighbor'.format(compare_with[i])
    
    # pdb.set_trace()
    
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    diff_array = []
    strct_from_array = []
    strct_to_array = []
    
    pvalue_array = []
    
    for s in range(len(type_strct)):
        
        xl_norm = pd.read_excel('{}/{}.xlsx'.format(xl_path,type_strct[s]),'Normal')
        
        cluster_type_list = np.array(xl_norm['Cluster Type'])
        
        xl_comp = pd.read_excel('{}/{}.xlsx'.format(xl_path,type_strct[s]),compare_with[i])
        
        norm_median = xl_norm['median']
        
        comp_median = xl_comp['median']
        
        diff_median = norm_median-comp_median
        
        if s == 0:
            diff_array = diff_median.copy()
        else:
            diff_array = pd.concat([diff_array,diff_median])
        
        for ss in range(len(xl_norm['Cluster Type'])):
            strct_from_array.append(xl_norm['Cluster Type'][ss])
            strct_to_array.append(type_strct[s].capitalize())
            cluster_idx = np.where(np.array(xl_norm['Cluster Type']) == xl_norm['Cluster Type'][ss])[0][0]
            res_stat = ranksums(xl_norm.iloc[cluster_idx][1:-1],xl_comp.iloc[cluster_idx][1:-1])
            pvalue_array.append(res_stat.pvalue)
        
        # pdb.set_trace()
    # pdb.set_trace()
    diff_array = np.array(diff_array).reshape((len(diff_array),1))
    strct_from_array = np.array(strct_from_array).reshape((len(strct_from_array),1))
    strct_to_array = np.array(strct_to_array).reshape((len(strct_to_array),1))
    pvalue_array = np.array(pvalue_array).reshape((len(pvalue_array),1))
    # pdb.set_trace()
    df_p = np.concatenate((strct_from_array,pvalue_array,strct_to_array),axis=1)
    
    df_p = pd.DataFrame(df_p)
    
    df_p.columns = ['Cluster Type', 'P-Value', 'Structure']
    
    df_pvalue = np.empty((len(cluster_type_list),len(type_strct)))
    df_pvalue = pd.DataFrame(df_pvalue)
    df_pvalue = df_pvalue.astype(str)
    
    # pdb.set_trace()
    
    df = np.concatenate((strct_from_array,diff_array,strct_to_array),axis=1)
    
    df = pd.DataFrame(df)
    
    df_forsave = np.concatenate((strct_from_array,strct_to_array,diff_array,pvalue_array),axis=1)
    
    df_forsave = pd.DataFrame(df_forsave)
    
    df_forsave.columns = ['Cluster Type','Structure', 'Delta', 'P-Value']
    
    df.columns = ['Cluster Type', 'Delta', 'Structure']
    
    df['Delta'] = df['Delta'].astype(float)
    
    df_forsave.to_excel('{}/Delta_Normal_{}.xlsx'.format(save_path, compare_with[i]), index=False)
    
    pivoted_df = df.pivot(index="Cluster Type", columns="Structure", values="Delta")
    
    for ii in range(len(pivoted_df.columns)):
        # pdb.set_trace()
        for ss in range(len(cluster_type_list)):
            p_idx = np.where((np.array(df_p['Structure']) == pivoted_df.columns[ii])&(np.array(df_p['Cluster Type']) == str(cluster_type_list[ss])))[0][0]
            
            piv_idx = np.where(np.array(pivoted_df.index) == str(cluster_type_list[ss]))[0][0]
            
            if float(df_p['P-Value'][p_idx]) > 0.05:
                df_pvalue.loc[piv_idx,ii] = ''
            else:
                df_pvalue.loc[piv_idx,ii] = '*'
    
    # pdb.set_trace()
    
    plt.figure(figsize=(30,30))
    sns.heatmap(pivoted_df,cmap='coolwarm',annot=df_pvalue,fmt='s',annot_kws={"size":50}, vmin=-1, vmax=1,center=0)
    plt.savefig("{}/Cellneighbor_heatmap_diff_N_{}.png".format(save_path,compare_with[i]))
    plt.close()
