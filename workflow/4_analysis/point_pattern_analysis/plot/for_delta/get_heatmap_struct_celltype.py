import os
import sys
import warnings
from glob import glob
import numpy as np
import pandas as pd
from scipy.stats import ranksums
import matplotlib.pyplot as plt
import seaborn as sns

# ---------- CONFIG ----------
input_dir = "path/to/data"
out_dir = "path/to/output"
MPN_SUBTYPES = ["Normal", "ET", "PV", "MF"]
STRUCTURES = ["arteriole", "bone", "fat", "sinusoid"]  # filename matching (case-insensitive)
IGNORED_SHEETS = {"All"}
os.makedirs(out_dir, exist_ok=True)

# ---------- HELPERS ----------
def sanitize_label(lbl):
    if pd.isna(lbl):
        return ""
    s = str(lbl).strip().replace("HSPC", "HSC")
    return " ".join(s.split())

def coerce_numeric_columns(df):
    df = df.copy()
    for col in df.columns:
        if col == "Cell Type":
            continue
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def find_rows_for_cell(df, cell_label):
    if df is None or df.empty:
        return None
    if "Cell Type" in df.columns:
        matches = df[df["Cell Type"].astype(str) == cell_label]
        if not matches.empty:
            return matches
    return None

def safe_vmax(arr):
    with np.errstate(all="ignore"):
        v = np.nanmax(np.abs(arr))
    return 1.0 if not np.isfinite(v) or v == 0 else float(v)

def build_cell_vs_cell_median(df_trans):
    """
    df_trans: dataframe with 'Cell Type' column and many numeric target columns.
    Returns: DataFrame index=source cell types, columns=target column names (sanitized), values=median
    """
    if "Cell Type" not in df_trans.columns:
        raise ValueError("Input df must contain 'Cell Type'")
    df_trans = df_trans.copy()
    df_trans["Cell Type"] = df_trans["Cell Type"].astype(str).apply(sanitize_label)
    src_labels = [x for x in df_trans["Cell Type"].unique() if x]
    candidate_cols = [c for c in df_trans.columns if c != "Cell Type"]
    display_cols = [sanitize_label(c) or str(c) for c in candidate_cols]
    mat = pd.DataFrame(index=src_labels, columns=display_cols, dtype=float)
    for src in src_labels:
        matches = find_rows_for_cell(df_trans, src)
        if matches is None or matches.empty:
            mat.loc[src, :] = np.nan
            continue
        for orig_col, disp_col in zip(candidate_cols, display_cols):
            vals = pd.to_numeric(matches[orig_col], errors="coerce").dropna().values
            mat.at[src, disp_col] = np.median(vals) if vals.size > 0 else np.nan
    return mat

# ---------- LOAD FILES & ASSIGN TO STRUCTURE ----------
print("Scanning input directory:", input_dir)
xlsx_paths = sorted(glob(os.path.join(input_dir, "*.xlsx")))
if not xlsx_paths:
    sys.exit("No Excel files found in the input directory.")

# per_structure_tables[structure][subtype] -> list of dataframes (one per file)
per_structure_tables = {s: {sub: [] for sub in MPN_SUBTYPES} for s in STRUCTURES}

for p in xlsx_paths:
    fname = os.path.basename(p).lower()
    matched_struct = None
    for s in STRUCTURES:
        if s in fname:
            matched_struct = s
            break
    if matched_struct is None:
        warnings.warn(f"File {p} did not match any known structure keywords; skipping.")
        continue
    try:
        xl = pd.ExcelFile(p)
    except Exception as e:
        warnings.warn(f"Cannot open {p}: {e}")
        continue
    for sub in MPN_SUBTYPES:
        if sub in xl.sheet_names and sub not in IGNORED_SHEETS:
            df = xl.parse(sub)
            if "Cell Type" not in df.columns:
                df = df.reset_index().rename(columns={"index": "Cell Type"})
            df["Cell Type"] = df["Cell Type"].astype(str).apply(sanitize_label)
            df = coerce_numeric_columns(df)
            per_structure_tables[matched_struct][sub].append(df)

# Ensure Normal exists for each structure (we allow missing but will warn)
for s in STRUCTURES:
    if not per_structure_tables[s]["Normal"]:
        warnings.warn(f"No Normal sheets found for structure '{s}'. Deltas for that structure will be NaN.")

# Concatenate lists into single DataFrame per structure/subtype
for s in STRUCTURES:
    for sub in MPN_SUBTYPES:
        if per_structure_tables[s][sub]:
            per_structure_tables[s][sub] = pd.concat(per_structure_tables[s][sub], ignore_index=True, sort=False)
        else:
            per_structure_tables[s][sub] = pd.DataFrame(columns=["Cell Type"])

# ---------- COMPUTE MEDIANS PER STRUCTURE × SUBTYPE ----------
median_per_structure = {s: {} for s in STRUCTURES}
for s in STRUCTURES:
    for sub in MPN_SUBTYPES:
        df = per_structure_tables[s][sub]
        if df.empty:
            median_per_structure[s][sub] = pd.DataFrame()
            continue
        mat = build_cell_vs_cell_median(df)
        # Save sorted (alphabetical) median tables for record
        mat_sorted = mat.reindex(index=sorted(mat.index, key=str.lower),
                                 columns=sorted(mat.columns, key=str.lower))
        median_per_structure[s][sub] = mat
        mat_sorted.to_excel(os.path.join(out_dir, f"Cell_vs_Cell_Median_{s}_{sub}.xlsx"))
        print(f"Saved median table for structure={s}, subtype={sub}")

# ---------- BUILD DELTA MATRICES (cols=structures) FOR ET/PV/MF vs Normal ----------
COMPARISONS = ["ET", "PV", "MF"]
# Collect union of all cell types and union of all target columns across all structures/subtypes
all_source_cells = set()
all_target_cols = set()
for s in STRUCTURES:
    for sub in MPN_SUBTYPES:
        mat = median_per_structure[s][sub]
        if not mat.empty:
            all_source_cells.update(mat.index.tolist())
            all_target_cols.update(mat.columns.tolist())

# But the user's requested heatmap: rows = source cell type, columns = structure.
# So we will compute delta per source cell type and per structure; target cell types are ignored for this axis.
# We will compute delta using the median value of the "self->target" entries aggregated by taking mean across target columns?
# Interpretation: If your median matrices are "source (rows) vs target (columns)", the user asked "heatmap should show cell type in y-axis and x-axis as each structure".
# Likely they want the aggregate (row-wise) median per source cell for that structure. We'll compute **row-wise median across target columns** for each source cell.
#
# Implementation: For each structure and subtype, compute row_median = median across target columns (ignoring NaN).
# Then delta(structure) = row_median(subtype) - row_median(Normal).
#
# Significance: for each (source cell, structure), perform ranksums comparing the underlying raw values used to compute medians:
# - For a given structure s, source cell r, and target column set: extract all numeric values from rows in the original per_structure_tables[s][sub] where 'Cell Type' == r (across all relevant target columns),
# - compare Normal vs subtype arrays with ranksums.

# Compute row-wise medians per structure/subtype
row_median_per_structure = {s: {} for s in STRUCTURES}
for s in STRUCTURES:
    for sub in MPN_SUBTYPES:
        mat = median_per_structure[s][sub]
        if mat.empty:
            row_median_per_structure[s][sub] = pd.Series(dtype=float)
            continue
        # row-wise median across all target columns
        row_med = mat.median(axis=1, skipna=True)
        # ensure index is sanitized strings
        row_med.index = [sanitize_label(i) for i in row_med.index]
        row_median_per_structure[s][sub] = row_med

# Build union of source cell types across all row_median series
all_sources = sorted({src for s in STRUCTURES for sub in MPN_SUBTYPES for src in row_median_per_structure[s][sub].index})

# Prepare results containers and compute significance
for compare_sub in COMPARISONS:
    # delta_df: rows = source cells, cols = structures
    delta_df = pd.DataFrame(index=all_sources, columns=[s.capitalize() for s in STRUCTURES], dtype=float)
    sig_df = pd.DataFrame("", index=all_sources, columns=[s.capitalize() for s in STRUCTURES], dtype=object)

    for s in STRUCTURES:
        normal_series = row_median_per_structure[s].get("Normal", pd.Series(dtype=float))
        comp_series = row_median_per_structure[s].get(compare_sub, pd.Series(dtype=float))

        for src in all_sources:
            val_norm = normal_series.get(src, np.nan)
            val_comp = comp_series.get(src, np.nan)
            if pd.isna(val_norm) and pd.isna(val_comp):
                delta = np.nan
            else:
                # compute delta (comp - norm) — user asked difference between Normal vs ET/PV/MF.
                # Interpreting as: delta = subtype - Normal
                delta = (val_comp if not pd.isna(val_comp) else np.nan) - (val_norm if not pd.isna(val_norm) else np.nan)
            delta_df.at[src, s.capitalize()] = delta

            # --- significance ---
            # Extract raw numeric values across all target columns for this source cell for Normal and compare_sub within this structure
            def extract_raw_for_source(struct_df, source_label):
                # struct_df is the concatenated dataframe for that structure and subtype
                if struct_df is None or struct_df.empty:
                    return np.array([])
                m = find_rows_for_cell(struct_df, source_label)
                if m is None or m.empty:
                    return np.array([])
                # collect all numeric values across all columns except 'Cell Type'
                vals = []
                for c in m.columns:
                    if c == "Cell Type":
                        continue
                    arr = pd.to_numeric(m[c], errors="coerce").dropna().tolist()
                    vals.extend(arr)
                return np.array(vals)

            arr_norm = extract_raw_for_source(per_structure_tables[s]["Normal"], src)
            arr_comp = extract_raw_for_source(per_structure_tables[s][compare_sub], src)

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

    # Save delta table and significance table
    delta_df.to_excel(os.path.join(out_dir, f"Delta_rowMedian_{compare_sub}_vs_Normal_by_structure.xlsx"))
    sig_df.to_excel(os.path.join(out_dir, f"Delta_rowMedian_signif_{compare_sub}_vs_Normal_by_structure.xlsx"))
    print(f"Saved delta and significance tables for {compare_sub} vs Normal")

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
        cmap="coolwarm", center=0, vmin=vmin, vmax=vmax,
        mask=mask, annot=sig_df.values, fmt="",
        annot_kws={"size": ANNOT_FONTSIZE, "weight": "bold"},
        cbar_kws={"shrink": 0.8}, ax=ax
    )

    ax.set_title(f"Δ {compare_sub} vs Normal — row-wise median per structure", fontsize=TITLE_FONTSIZE, pad=30)
    ax.set_xlabel("Structure", fontsize=LABEL_FONTSIZE, labelpad=20)
    ax.set_ylabel("Source Cell Type", fontsize=LABEL_FONTSIZE, labelpad=20)
    ax.tick_params(axis='x', labelsize=TICK_FONTSIZE, rotation=0)
    ax.tick_params(axis='y', labelsize=TICK_FONTSIZE)
    # colorbar font
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=COLORBAR_FONTSIZE)

    plt.tight_layout()
    out_png = os.path.join(out_dir, f"Delta_rowMedian_{compare_sub}_vs_Normal_by_structure.png")
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved heatmap: {out_png}")

print("DONE. All outputs saved in:", os.path.abspath(out_dir))
