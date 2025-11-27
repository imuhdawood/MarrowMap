type_strct = ['Bone','Fat','Arteriole','Sinusoid']

xl_path = 'xenium_pipeline/ppm_strct_strct_patch_more'

# xl_path = 'xenium_pipeline/gathered_ppm_downsampled_cell/Granulocyte_mast.xlsx'

compare_with = ['MF','ET','PV']


for i in range(len(compare_with)):
    save_path = 'xenium_pipeline/delta_norm{}_strct_to_strct'.format(compare_with[i])

    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    diff_array = []
    strct_from_array = []
    strct_to_array = []
    
    pvalue_array = []
    
    for s in range(len(type_strct)):
        
        xl_norm = pd.read_excel('{}/Normal.xlsx'.format(xl_path),type_strct[s])
        
        for ss in range(len(type_strct)):
            
            strct_idx = np.where(np.array(xl_norm['structure']) == type_strct[ss])[0][0]
            
            # pdb.set_trace()
            
            xl_comp = pd.read_excel('{}/{}.xlsx'.format(xl_path,compare_with[i]),type_strct[s])
            
            # pdb.set_trace()
            
            diff_array.append(xl_norm['median'][strct_idx]-xl_comp['median'][strct_idx])
            
            strct_from_array.append(type_strct[ss])
            
            strct_to_array.append(type_strct[s])
            
            res_stat = ranksums(xl_norm.iloc[strct_idx][1:-1],xl_comp.iloc[strct_idx][1:-1])
            # res_stat = ranksums(xl_n.iloc[cc][1:],xl_mf.iloc[mf_idx][1:],alternative='greater')
            pvalue_array.append(res_stat.pvalue)
            
            # pdb.set_trace()
        
        # pdb.set_trace()
    # pdb.set_trace()
    diff_array = np.array(diff_array).reshape((len(diff_array),1))
    strct_from_array = np.array(strct_from_array).reshape((len(strct_from_array),1))
    strct_to_array = np.array(strct_to_array).reshape((len(strct_to_array),1))
    pvalue_array = np.array(pvalue_array).reshape((len(pvalue_array),1))
    
    df_p = np.concatenate((strct_from_array,pvalue_array,strct_to_array),axis=1)
    
    df_p = pd.DataFrame(df_p)
    
    df_p.columns = ['Structure From', 'P-Value', 'Structure To']
    
    df_pvalue = np.empty((len(type_strct),len(type_strct)))
    df_pvalue = pd.DataFrame(df_pvalue)
    df_pvalue = df_pvalue.astype(str)
    
    df = np.concatenate((strct_from_array,diff_array,strct_to_array),axis=1)
    
    df = pd.DataFrame(df)
    
    df_forsave = np.concatenate((strct_from_array,strct_to_array,diff_array,pvalue_array),axis=1)
    
    df_forsave = pd.DataFrame(df_forsave)
    
    df_forsave.columns = ['Structure From','Structure To', 'Delta', 'P-Value']
    
    df.columns = ['Structure From', 'Delta', 'Structure To']
    
    df['Delta'] = df['Delta'].astype(float)
    
    df_forsave.to_excel('{}/Delta_Normal_{}.xlsx'.format(save_path, compare_with[i]), index=False)
    
    pivoted_df = df.pivot(index="Structure From", columns="Structure To", values="Delta")
    
    for ii in range(len(pivoted_df.columns)):
        # pdb.set_trace()
        for ss in range(len(type_strct)):
            p_idx = np.where((np.array(df_p['Structure To']) == pivoted_df.columns[ii])&(np.array(df_p['Structure From']) == type_strct[ss]))[0][0]
            
            piv_idx = np.where(np.array(pivoted_df.index) == type_strct[ss])[0][0]
            
            if float(df_p['P-Value'][p_idx]) > 0.05:
                df_pvalue.loc[piv_idx,ii] = ''
            else:
                df_pvalue.loc[piv_idx,ii] = '*'
    
    # pdb.set_trace()
    
    plt.figure(figsize=(30,30))
    sns.heatmap(pivoted_df,cmap='coolwarm',annot=df_pvalue,fmt='s',annot_kws={"size":50}, vmin=-1, vmax=1,center=0)
    plt.savefig("{}/Structure_heatmap_diff_N_{}.png".format(save_path,compare_with[i]))
    plt.close()
