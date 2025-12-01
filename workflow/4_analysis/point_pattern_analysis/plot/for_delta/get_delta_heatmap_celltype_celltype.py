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
    numeric_cols = []
    for col in df.columns:
        if col == "Cell Type":
            continue
        coerced = pd.to_numeric(df[col], errors="coerce")
        if coerced.notna().any():
            df[col] = coerced
            numeric_cols.append(col)
    return df, numeric_cols

def find_rows_for_cell(df, cell_label):
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
    if "Cell Type" not in df_trans.columns:
        raise ValueError("Input df must contain 'Cell Type'")
    df_trans = df_trans.copy()
    df_trans["Cell Type"] = df_trans["Cell Type"].astype(str).apply(sanitize_label)
    src_labels = [x for x in df_trans["Cell Type"].unique() if x]
    candidate_cols = []
    for c in df_trans.columns:
        if c == "Cell Type":
            continue
        disp = sanitize_label(c)
        coerced = pd.to_numeric(df_trans[c], errors="coerce")
        if coerced.notna().any() or (disp in src_labels):
            candidate_cols.append(c)
    display_cols = [sanitize_label(c) or str(c) for c in candidate_cols]
    mat = pd.DataFrame(index=src_labels, columns=display_cols, dtype=float)
    for src in src_labels:
        matches = find_rows_for_cell(df_trans, src)
        if matches is None or matches.empty:
            mat.loc[src, :] = np.nan
            continue
        for orig_col, disp_col in zip(candidate_cols, display_cols):
            vals = pd.to_numeric(matches[orig_col], errors="coerce").values
            nums = [float(v) for v in vals if pd.notna(v) and np.isfinite(v)]
            mat.at[src, disp_col] = np.median(nums) if nums else np.nan
    return mat

# ---------- LOAD ----------
print("Scanning input directory:", input_dir)
xlsx_paths = sorted(glob(os.path.join(input_dir, "*.xlsx")))
if not xlsx_paths:
    sys.exit("No Excel files found in the input directory.")

per_subtype_tables = {s: [] for s in MPN_SUBTYPES}

for p in xlsx_paths:
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
            df, _ = coerce_numeric_columns(df)
            per_subtype_tables[sub].append(df)

if not per_subtype_tables["Normal"]:
    sys.exit("No Normal sheets found; cannot compute delta.")

for sub in MPN_SUBTYPES:
    if per_subtype_tables[sub]:
        per_subtype_tables[sub] = pd.concat(per_subtype_tables[sub], ignore_index=True, sort=False)
    else:
        per_subtype_tables[sub] = pd.DataFrame(columns=["Cell Type"])

# ---------- MEDIAN MATRICES ----------
median_matrices = {}
for sub, df in per_subtype_tables.items():
    if df.empty:
        median_matrices[sub] = pd.DataFrame()
        continue
    mat = build_cell_vs_cell_median(df)

    # Sort alphabetically only for final saved output
    mat_sorted = mat.copy()
    mat_sorted = mat_sorted.reindex(index=sorted(mat.index, key=str.lower),
                                    columns=sorted(mat.columns, key=str.lower))

    median_matrices[sub] = mat
    mat_sorted.to_excel(os.path.join(out_dir, f"Cell_vs_Cell_Median_{sub}.xlsx"))
    print("Saved median table for", sub)

# ---------- DELTA vs NORMAL ----------
def build_display_to_original_map(df):
    m = {}
    for c in df.columns:
        if c == "Cell Type":
            continue
        disp = sanitize_label(c) or str(c)
        m.setdefault(disp, []).append(c)
    return m

map_norm = build_display_to_original_map(per_subtype_tables["Normal"])
normal_mat = median_matrices["Normal"]

for sub in ["ET", "PV", "MF"]:
    df_sub = per_subtype_tables[sub]
    if df_sub.empty:
        print(f"Skipping {sub}, no data.")
        continue
    sub_mat = median_matrices[sub]
    rows = list(set(normal_mat.index) | set(sub_mat.index))
    cols = list(set(normal_mat.columns) | set(sub_mat.columns))
    norm_aligned = pd.DataFrame(index=rows, columns=cols)
    sub_aligned = pd.DataFrame(index=rows, columns=cols)
    for r in rows:
        for c in cols:
            norm_aligned.at[r, c] = normal_mat.at[r, c] if (r in normal_mat.index and c in normal_mat.columns) else np.nan
            sub_aligned.at[r, c] = sub_mat.at[r, c] if (r in sub_mat.index and c in sub_mat.columns) else np.nan
    delta = sub_aligned - norm_aligned
    sig = pd.DataFrame("", index=rows, columns=cols)
    map_sub = build_display_to_original_map(df_sub)

    def extract_values(tbl, mapping, src, tgt):
        m = find_rows_for_cell(tbl, src)
        if m is None or m.empty:
            return np.array([])
        cols = mapping.get(tgt, [tgt])
        vals = []
        for c in cols:
            if c not in m.columns:
                continue
            arr = pd.to_numeric(m[c], errors="coerce").dropna().to_list()
            vals.extend(arr)
        return np.array(vals)

    for r in rows:
        for c in cols:
            arr_n = extract_values(per_subtype_tables["Normal"], map_norm, r, c)
            arr_s = extract_values(df_sub, map_sub, r, c)
            if len(arr_n) and len(arr_s) and (np.std(arr_n) > 0 or np.std(arr_s) > 0):
                try:
                    p = ranksums(arr_n, arr_s).pvalue
                    if p <= 0.05:
                        sig.at[r, c] = "*"
                except Exception:
                    pass

    # Sort alphabetically only when saving
    delta_sorted = delta.reindex(index=sorted(rows, key=str.lower),
                                 columns=sorted(cols, key=str.lower))
    sig_sorted = sig.reindex(index=sorted(rows, key=str.lower),
                             columns=sorted(cols, key=str.lower))

    delta_sorted.to_excel(os.path.join(out_dir, f"Delta_{sub}_vs_Normal_Cell_vs_Cell.xlsx"))
    print(f"Saved delta table for {sub} vs Normal")

    # Plot delta heatmap
    vals = delta_sorted.astype(float).values
    mask = np.isnan(vals)
    if np.all(mask):
        continue
    vmax = safe_vmax(vals)
    annot = np.where(mask, "", sig_sorted.values)
    fig, ax = plt.subplots(figsize=(max(6, vals.shape[1]*0.35), max(4, vals.shape[0]*0.25)))
    sns.heatmap(pd.DataFrame(vals, index=delta_sorted.index, columns=delta_sorted.columns),
                cmap="coolwarm", center=0, vmin=-vmax, vmax=vmax,
                mask=mask, annot=annot, fmt="", annot_kws={"size": 10}, ax=ax)
    ax.set(title=f"Delta {sub} vs Normal — Cell vs Cell",
           xlabel="Target Cell Type", ylabel="Source Cell Type")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"Delta_heatmap_{sub}_vs_Normal.png"), dpi=150)
    plt.close(fig)

print("DONE. Outputs saved in:", os.path.abspath(out_dir))
