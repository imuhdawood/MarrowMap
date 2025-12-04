import os
import sys
import warnings
from glob import glob
import numpy as np
import pandas as pd
from scipy.stats import ranksums
import matplotlib.pyplot as plt
import seaborn as sns

compare_type = ['struct_struct', 'struct_celltype', 'struct_cellneighbor',
                'celltype_celltype','cellneighbor_cellneighbor']

MPN_SUBTYPES = ["Normal", "ET", "PV", "MF"]

COMPARISONS = ["ET", "PV", "MF"]

def sanitize_label(lbl):
    if pd.isna(lbl):
        return ""
    s = str(lbl).strip().replace("HSPC", "HSC")
    return " ".join(s.split())

def coerce_numeric_columns(df):
    df = df.copy()
    # 
    for col in df.columns:
        if col == df.columns[0]:
            continue
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def find_rows_for_from(df, src_labels):

    if df is None or df.empty:
        return None
    
    matches = df[df[df.columns[0]].astype(str) == src_labels]
    if not matches.empty:
        return matches
    
    
    return None

def safe_vmax(arr):
    with np.errstate(all="ignore"):
        v = np.nanmax(np.abs(arr))
    return 1.0 if not np.isfinite(v) or v == 0 else float(v)

def build_median(df_trans):
    """
    df_trans: dataframe with 'Cell Type' column and many numeric target columns.
    Returns: DataFrame index=source cell types, columns=target column names (sanitized), values=median
    """
    
    df_trans = df_trans.copy()
    df_trans[df_trans.columns[0]] = df_trans[df_trans.columns[0]].astype(str).apply(sanitize_label)
    src_labels = [x for x in df_trans[df_trans.columns[0]].unique() if x]
    candidate_cols = [c for c in df_trans.columns if c != df_trans.columns[0]]
    display_cols = [sanitize_label(c) or str(c) for c in candidate_cols]
    mat = pd.DataFrame(index=src_labels, columns=display_cols, dtype=float)
    # 
    for src in src_labels:
        matches = find_rows_for_from(df_trans, src)
        if matches is None or matches.empty:
            mat.loc[src, :] = np.nan
            continue
        for orig_col, disp_col in zip(candidate_cols, display_cols):
            vals = pd.to_numeric(matches[orig_col], errors="coerce").dropna().values
            mat.at[src, disp_col] = np.median(vals) if vals.size > 0 else np.nan
    return mat

for comptype in compare_type:
    
    if 'struct' in comptype:
        To_grp = ["Arteriole", "Bone", "Fat", "Sinusoid"]
        if comptype == 'struct_struct':
            # ---------- CONFIG ----------
            input_dir = "../gather_data/{}".format(comptype)
            out_dir = "../delta/{}".format(comptype)
            IGNORED_SHEETS = {"All"}
            os.makedirs(out_dir, exist_ok=True)
            
            # ---------- LOAD FILES & ASSIGN TO STRUCTURE ----------
            print("Scanning input directory:", input_dir)
            xlsx_paths = sorted(glob(os.path.join(input_dir, "*.xlsx")))
            if not xlsx_paths:
                sys.exit("No Excel files found in the input directory.")
            
            per_to_tables = {s: {sub: [] for sub in MPN_SUBTYPES} for s in To_grp}
            # 
            for sub in MPN_SUBTYPES:
                xl = pd.ExcelFile('{}/{}.xlsx'.format(input_dir,sub))
                
                for tg in To_grp:
                    if tg in xl.sheet_names:
                        df = xl.parse(tg)
                        
                        df[df.columns[0]] = df[df.columns[0]].astype(str)
                        df = coerce_numeric_columns(df)
                        per_to_tables[tg][sub].append(df)
            
            
            
        elif (comptype == 'struct_celltype') | (comptype == 'struct_cellneighbor'):
    
            # ---------- CONFIG ----------
            input_dir = "../gather_data/{}_".format(comptype)
            out_dir = "../delta/{}_".format(comptype)
            IGNORED_SHEETS = {"All"}
            os.makedirs(out_dir, exist_ok=True)
            
            
            # ---------- LOAD FILES & ASSIGN TO STRUCTURE ----------
            print("Scanning input directory:", input_dir)
            xlsx_paths = sorted(glob(os.path.join(input_dir, "*.xlsx")))
            if not xlsx_paths:
                sys.exit("No Excel files found in the input directory.")
            
            per_to_tables = {s: {sub: [] for sub in MPN_SUBTYPES} for s in To_grp}
            # 
            for p in xlsx_paths:
                fname = os.path.basename(p)
                matched_togrp = None
                for s in To_grp:
                    if s in fname:
                        matched_togrp = s
                        break
                if matched_togrp is None:
                    warnings.warn(f"File {p} did not match any known structure keywords; skipping.")
                    continue
                try:
                    xl = pd.ExcelFile(p)
                except Exception as e:
                    warnings.warn(f"Cannot open {p}: {e}")
                    continue
                for sub in MPN_SUBTYPES:
                    # 
                    if sub in xl.sheet_names and sub not in IGNORED_SHEETS:
                        df = xl.parse(sub)
                        
                        df[df.columns[0]] = df[df.columns[0]].astype(str).apply(sanitize_label)
                        df = coerce_numeric_columns(df)
                        per_to_tables[matched_togrp][sub].append(df)
            
            # Ensure Normal exists for each structure (we allow missing but will warn)
            for s in To_grp:
                if not per_to_tables[s]["Normal"]:
                    warnings.warn(f"No Normal sheets found for structure '{s}'. Deltas for that structure will be NaN.")
            
        
        else:
            raise ValueError("For Strucrure, only accepts [struct, celltype, cellneighbor].")
        
        # 
        
    else:
        
        if comptype == 'celltype_celltype':
            To_grp = ['Adipo-MSC','B_cell', 'DC', 'Endothelial', 'Erythroid', 'GMP','Granulocyte/mast', 'HSPC',
                          'Macrophage', 'Megakaryocyte', 'Monocyte', 'Myeloid', 'Osteo-MSC', 'Plasma_cell', 
                          'SMC', 'Stromal', 'T_cell']
            
            # ---------- CONFIG ----------
            input_dir = "../gather_data/{}".format(comptype)
            out_dir = "../delta/{}".format(comptype)
            IGNORED_SHEETS = {"All"}
            os.makedirs(out_dir, exist_ok=True)
            
            # ---------- LOAD FILES & ASSIGN TO STRUCTURE ----------
            print("Scanning input directory:", input_dir)
            xlsx_paths = sorted(glob(os.path.join(input_dir, "*.xlsx")))
            if not xlsx_paths:
                sys.exit("No Excel files found in the input directory.")
            
            per_to_tables = {s: {sub: [] for sub in MPN_SUBTYPES} for s in To_grp}
            
            for p in xlsx_paths:
                fname = os.path.basename(p)
                if fname.split('.xlsx')[0] == 'Granulocyte_mast':
                    fname = '{}/{}.xlsx'.format(fname.split('_')[0],fname.split('.')[0].split('_')[1])
                matched_togrp = None
                
                for s in To_grp:
                    if s in fname:
                        matched_togrp = s
                        break
                if matched_togrp is None:
                    warnings.warn(f"File {p} did not match any known structure keywords; skipping.")
                    continue
                try:
                    xl = pd.ExcelFile(p)
                except Exception as e:
                    warnings.warn(f"Cannot open {p}: {e}")
                    continue
                for sub in MPN_SUBTYPES:
                    # 
                    if sub in xl.sheet_names and sub not in IGNORED_SHEETS:
                        df = xl.parse(sub)
                        
                        df[df.columns[0]] = df[df.columns[0]].astype(str).apply(sanitize_label)
                        df = coerce_numeric_columns(df)
                        per_to_tables[matched_togrp][sub].append(df)
                        
        elif comptype == 'cellneighbor_cellneighbor':
            To_grp = ['0','1','2','3','4','5','6','7','8','9']
            
            # ---------- CONFIG ----------
            input_dir = "../gather_data/{}".format(comptype)
            out_dir = "../delta/{}".format(comptype)
            IGNORED_SHEETS = {"All"}
            os.makedirs(out_dir, exist_ok=True)
            
            # ---------- LOAD FILES & ASSIGN TO STRUCTURE ----------
            print("Scanning input directory:", input_dir)
            xlsx_paths = sorted(glob(os.path.join(input_dir, "*.xlsx")))
            if not xlsx_paths:
                sys.exit("No Excel files found in the input directory.")
            
            per_to_tables = {s: {sub: [] for sub in MPN_SUBTYPES} for s in To_grp}
            # 
            for p in xlsx_paths:
                fname = os.path.basename(p)
                matched_togrp = None
                # 
                for s in To_grp:
                    if s in fname:
                        matched_togrp = s
                        break
                if matched_togrp is None:
                    warnings.warn(f"File {p} did not match any known structure keywords; skipping.")
                    continue
                try:
                    xl = pd.ExcelFile(p)
                except Exception as e:
                    warnings.warn(f"Cannot open {p}: {e}")
                    continue
                for sub in MPN_SUBTYPES:
                    # 
                    if sub in xl.sheet_names and sub not in IGNORED_SHEETS:
                        df = xl.parse(sub)
                        
                        df[df.columns[0]] = df[df.columns[0]].astype(str).apply(sanitize_label)
                        df = coerce_numeric_columns(df)
                        per_to_tables[matched_togrp][sub].append(df)
            
    # 
    # Concatenate lists into single DataFrame per structure/subtype
    for s in To_grp:
        for sub in MPN_SUBTYPES:
            if per_to_tables[s][sub]:
                per_to_tables[s][sub] = pd.concat(per_to_tables[s][sub], ignore_index=True, sort=False)
            else:
                per_to_tables[s][sub] = pd.DataFrame(columns=[df.columns[0]])
    
    # ---------- COMPUTE MEDIANS PER STRUCTURE × SUBTYPE ----------
    median_per_togrp = {s: {} for s in To_grp}
    for s in To_grp:
        for sub in MPN_SUBTYPES:
            df = per_to_tables[s][sub]
            if df.empty:
                median_per_togrp[s][sub] = pd.DataFrame()
                continue
            
            mat = build_median(df)
            # Save sorted (alphabetical) median tables for record
            mat_sorted = mat.reindex(index=sorted(mat.index, key=str.lower),
                                     columns=sorted(mat.columns, key=str.lower))
            median_per_togrp[s][sub] = mat
            print(f"Saved median table for structure={s}, subtype={sub}")
    
    row_median_per_togrp = {s: {} for s in To_grp}
    for s in To_grp:
        for sub in MPN_SUBTYPES:
            mat = median_per_togrp[s][sub]
            if mat.empty:
                row_median_per_togrp[s][sub] = pd.Series(dtype=float)
                continue
            # row-wise median across all target columns
            row_med = mat.median(axis=1, skipna=True)
            # ensure index is sanitized strings
            row_med.index = [sanitize_label(i) for i in row_med.index]
            row_median_per_togrp[s][sub] = row_med
    
    # Build union of source cell types across all row_median series
    all_sources = sorted({src for s in To_grp for sub in MPN_SUBTYPES for src in row_median_per_togrp[s][sub].index})
    
    # Prepare results containers and compute significance
    for compare_sub in COMPARISONS:
        # delta_df: rows = source cells, cols = structures
        delta_df = pd.DataFrame(index=all_sources, columns=[s.capitalize() for s in To_grp], dtype=float)
        sig_df = pd.DataFrame("", index=all_sources, columns=[s.capitalize() for s in To_grp], dtype=object)
    
        for s in To_grp:
            normal_series = row_median_per_togrp[s].get("Normal", pd.Series(dtype=float))
            comp_series = row_median_per_togrp[s].get(compare_sub, pd.Series(dtype=float))
    
            for src in all_sources:
                val_norm = normal_series.get(src, np.nan)
                val_comp = comp_series.get(src, np.nan)
                # 
                if pd.isna(val_norm) and pd.isna(val_comp):
                    delta = np.nan
                else:
                    # compute delta (comp - norm) — user asked difference between Normal vs ET/PV/MF.
                    # Interpreting as: delta = subtype - Normal
                    delta = (val_norm if not pd.isna(val_norm) else np.nan) - (val_comp if not pd.isna(val_comp) else np.nan)
                delta_df.at[src, s.capitalize()] = delta
    
                # --- significance ---
                # Extract raw numeric values across all target columns for this source cell for Normal and compare_sub within this structure
                def extract_raw_for_source(df_, source_label):
                    # struct_df is the concatenated dataframe for that structure and subtype
                    if df_ is None or df_.empty:
                        return np.array([])
                    # 
                    m = find_rows_for_from(df_, source_label)
                    if m is None or m.empty:
                        return np.array([])
                    # collect all numeric values across all columns except 'Cell Type'
                    vals = []
                    for c in m.columns:
                        if c == m.columns[0]:
                            continue
                        arr = pd.to_numeric(m[c], errors="coerce").dropna().tolist()
                        vals.extend(arr)
                    return np.array(vals)
    
                arr_norm = extract_raw_for_source(per_to_tables[s]["Normal"], src)
                arr_comp = extract_raw_for_source(per_to_tables[s][compare_sub], src)
                
                arr_norm = arr_norm[:-1]
                arr_comp = arr_comp[:-1]
                
                pval = None
                try:
                    if arr_norm.size and arr_comp.size and (np.std(arr_norm) > 0 or np.std(arr_comp) > 0):
                        pval = ranksums(arr_norm, arr_comp).pvalue
                    else:
                        pval = np.nan
                except Exception:
                    pval = np.nan
    
                if pval is not None and not np.isnan(pval) and pval <= 0.05:
                    sig_df.at[src, s.capitalize()] = "*"
                else:
                    sig_df.at[src, s.capitalize()] = ""
    
        # ---------- PLOT HEATMAP ----------
        vals = delta_df.astype(float).values
        mask = np.isnan(vals)
        if np.all(mask):
            print(f"No data for {compare_sub} vs Normal (all NaN). Skipping plot.")
            continue
    
        # Font scaling for figsize 30x30
        FIGSIZE = (30, 30)
        TITLE_FONTSIZE = 46
        LABEL_FONTSIZE = 40
        TICK_FONTSIZE = 30
        ANNOT_FONTSIZE = 36
        COLORBAR_FONTSIZE = 28
    
        vmax = safe_vmax(vals)
        # center at 0
        vmin = -vmax
        vmax = vmax
    
        fig, ax = plt.subplots(figsize=FIGSIZE)
        sns.heatmap(
            pd.DataFrame(vals, index=delta_df.index, columns=delta_df.columns),
            cmap="coolwarm", center=0, vmin=-1, vmax=1,
            mask=mask, annot=sig_df.values, fmt="",
            annot_kws={"size": ANNOT_FONTSIZE, "weight": "bold"},
            cbar_kws={"shrink": 0.8}, ax=ax
        )
    
        ax.set_title(f"Δ {compare_sub} vs Normal — row-wise median per structure", fontsize=TITLE_FONTSIZE, pad=30)
        ax.set_xlabel("Structure", fontsize=LABEL_FONTSIZE, labelpad=20)
        ax.set_ylabel("Source {}".format(df.columns[0]), fontsize=LABEL_FONTSIZE, labelpad=20)
        ax.tick_params(axis='x', labelsize=TICK_FONTSIZE, rotation=0)
        ax.tick_params(axis='y', labelsize=TICK_FONTSIZE)
        # colorbar font
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=COLORBAR_FONTSIZE)
    
        plt.tight_layout()
        out_png = os.path.join(out_dir, f"Delta_{compare_sub}_vs_Normal_by_structure.png")
        plt.savefig(out_png, dpi=300)
        plt.close(fig)
        print(f"Saved heatmap: {out_png}")
    
    print("DONE. All outputs saved in:", os.path.abspath(out_dir))
