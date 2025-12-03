import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import statsmodels.stats.multitest as smm
import seaborn as sns
import matplotlib.pyplot as plt
import os


def gather_interaction_values(
    mDf, 
    df, 
    dkey='Diagnosis_2', 
    group_labels=('Normal','MF'),
    ITS=False
):
    """
    Gathers per-sample interaction values for each cell-type pair, 
    separating them into the two groups (Normal vs MF).
    
    Parameters
    ----------
    mDf : pd.DataFrame
        Metadata DataFrame with a column dkey indicating the group 
        each sample belongs to (e.g., 'Normal' or 'MF').
    df : pd.DataFrame
        DataFrame that contains NxN interaction info for each sample. 
        We assume it has a column 'index' identifying the sample 
        (like '10693-Trefine', '10693-ITS', etc.).
    dkey : str
        Column in mDf that tells us which group the sample belongs to.
    group_labels : tuple
        A tuple or list with the two group labels we are interested in (e.g. ("Normal","MF")).
    ITS : bool
        Whether to select only ITS or Trefine rows from df. 
        (Adapts the subset logic from your original code.)

    Returns
    -------
    interaction_values : dict
        A nested dict of the form:
        {
          (cell_i, cell_j): {
              "Normal": [val1, val2, ...],
              "MF": [val1, val2, ...]
          },
          ...
        }
    cell_types : list
        The list of cell types used.
    """
    
    group_normal, group_mf = group_labels  # e.g. "Normal", "MF"
    
    # Identify which samples belong to each group
    normal_samples = set(mDf[mDf[dkey] == group_normal].index)
    mf_samples = set(mDf[mDf[dkey] == group_mf].index)
    
    # Subset df based on ITS / Trefine
    if ITS:
        df = df[~df['index'].str.contains('Trefine')]
    else:
        df = df[df['index'].str.contains('Trefine')]

    # Identify columns that correspond to cell types
    exclude_cols = ['index', 'CD69']
    cell_types = [c for c in df.columns if c not in exclude_cols]

    # Creating a dictionary to hold lists of interaction values
    # for each group, for each cell-type pair.
    interaction_values = {
        (i, j): {group_normal: [], group_mf: []} 
        for i in cell_types 
        for j in cell_types
    }
    
    # For each sample (ITS, or individual) in df, we find the NxN matrix and add values
    samples_in_df = df['index'].unique()
    
    for skey in samples_in_df:

        #Sample Selection
        short_key = skey.split('-')[0]
        
        if short_key in normal_samples:
            group_used = group_normal
        elif short_key in mf_samples:
            group_used = group_mf
        else:
            continue  # skip if not part of any of the group
        
        # Subset the DataFrame for this sample's to get CeLL X CELL Interaction matrix
        sample_block = df[df['index'] == skey].reindex(columns=cell_types, index=cell_types).fillna(0)
        
        # Now, accumulate each cell-type pair's value
        for i in cell_types:
            for j in cell_types:
                val_ij = sample_block.loc[i, j]
                interaction_values[(i, j)][group_used].append(val_ij)
    
    return interaction_values, cell_types

def compute_mannwhitney(
    interaction_values, 
    group_labels=('Normal','MF'),
    alpha=0.05
):
    """
    For each (i,j) pair, run a Mann-Whitney U test comparing the 
    distributions in the two groups. Correct for multiple tests via FDR (BH).
    
    Returns
    -------
    pval_df : pd.DataFrame
        DataFrame of BH-corrected p-values, indexed by (cell_i, cell_j).
    """
    
    group_normal, group_mf = group_labels
    
    pairs = sorted(interaction_values.keys())
    raw_pvals = []
    
    for pair in pairs:
        vals_normal = interaction_values[pair][group_normal]
        vals_mf = interaction_values[pair][group_mf]
        
        # If either group has < 2 data points, skip
        if len(vals_normal) < 2 or len(vals_mf) < 2:
            raw_pvals.append(np.nan)
            continue
        
        # Mann-Whitney U test (two-sided)
        stat, pval = mannwhitneyu(vals_normal, vals_mf, alternative='two-sided')
        raw_pvals.append(pval)
    
    raw_pvals = np.array(raw_pvals, dtype=float)
    
    # multiple test correction (FDR, Benjamini-Hochberg)
    mask = ~np.isnan(raw_pvals)
    rejected, corrected, _, _ = smm.multipletests(raw_pvals[mask], alpha=alpha, method='fdr_bh')
    
    # Place corrected pvals back into a results array
    corrected_pvals = np.full_like(raw_pvals, np.nan, dtype=float)
    corrected_pvals[mask] = corrected
    
    # Build a DataFrame of corrected p-values
    pval_df = pd.DataFrame(index=pairs, data={'pval_corrected': corrected_pvals})
    return pval_df

def compute_rank_biserial(vals_x, vals_y):
    """
    Compute rank-biserial correlation from Mann-Whitney U for two arrays.
    """
    # Mann-Whitney U
    U2, _ = mannwhitneyu(vals_x, vals_y, alternative='two-sided') # Scipy returns U1 related to X
    nx = len(vals_x)
    ny = len(vals_y)
    # rank-biserial correlation
    rb = 1.0 - (2.0 * U2) / (nx * ny)
    #rb = ((2.0 * U1) / (nx * ny))-1
    return rb

def compute_summary_matrices(interaction_values, cell_types, group_labels=('Normal','MF')):
    """
    Computes summary statistics (mean and median) for each cell-cell interaction
    matrix across samples in each group.
    
    Parameters
    ----------
    interaction_values : dict
        Nested dictionary with keys as (cell_i, cell_j) and values as dictionaries 
        containing lists of interactions for each group.
    cell_types : list
        List of cell types to use for the matrix axes.
    group_labels : tuple
        Tuple containing the two group labels (e.g., ("Normal", "MF")).
    
    Returns
    -------
    summary : dict
        A dictionary where each key is a group label and its value is another
        dict with two keys 'mean' and 'median' that contain the corresponding 
        DataFrames (cell_types x cell_types).
    """
    summary = {}
    for group in group_labels:
        mean_mat = pd.DataFrame(index=cell_types, columns=cell_types, dtype=float)
        median_mat = pd.DataFrame(index=cell_types, columns=cell_types, dtype=float)
        for i in cell_types:
            for j in cell_types:
                values = interaction_values[(i, j)][group]
                # Compute mean and median if values exist, else NaN
                if values:
                    mean_val = np.mean(values)
                    median_val = np.median(values)
                else:
                    mean_val = np.nan
                    median_val = np.nan
                mean_mat.loc[i, j] = mean_val
                median_mat.loc[i, j] = median_val
        summary[group] = {'mean': mean_mat, 'median': median_mat}
    return summary


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

EXP_ID = 'spatial_stats_RC_22122024_final_osteo'

FILE_PATH = f'/{EXP_ID}/interactions.csv' ## interaction matrix computed by running squidpy interaction matrix

META_FILE = '/xenium_keys.csv'

SUBTYPE1 = 'Normal'
SUBTYPE2 = 'PV'
ITS = False
dkey = 'Diagnosis_2'

OUTPUT_DIR = f'/OUTPUT/{EXP_ID}/Interactions_Plot'

def mkdir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


df = pd.read_csv(FILE_PATH, index_col='Unnamed: 0')
mDf = pd.read_csv(META_FILE, index_col="sample_key")

# ------------------------------------------------------------------------
# 1) Gather per-sample distributions
# ------------------------------------------------------------------------
interaction_values, cell_types = gather_interaction_values(
    mDf, df,
    dkey='Diagnosis_2',
    group_labels=(SUBTYPE1,SUBTYPE2),
    ITS=ITS
)

# ------------------------------------------------------------------------
# 2) Compute corrected p-values for each (i, j)
# ------------------------------------------------------------------------
pval_df = compute_mannwhitney(
    interaction_values,
    group_labels=(SUBTYPE1,SUBTYPE2),
    alpha=0.05
)

# ------------------------------------------------------------------------
# 3) Compute rank-biserial effect sizes
# ------------------------------------------------------------------------
effect_sizes = {}
for pair in pval_df.index:
    i, j = pair
    vals_normal = interaction_values[pair][SUBTYPE1]
    vals_mf     = interaction_values[pair][SUBTYPE2]
    if len(vals_normal) >= 2 and len(vals_mf) >= 2:
        effect_sizes[pair] = compute_rank_biserial(vals_normal, vals_mf)
    else:
        effect_sizes[pair] = np.nan

# ------------------------------------------------------------------------
# 4) Build DataFrames shaped (cell_types x cell_types) for heatmaps
# ------------------------------------------------------------------------
pval_mat = pd.DataFrame(index=cell_types, columns=cell_types, dtype=float)
effect_mat = pd.DataFrame(index=cell_types, columns=cell_types, dtype=float)

for (i, j), row in pval_df.iterrows():
    pval_mat.loc[i, j]   = row['pval_corrected']
    effect_mat.loc[i, j] = effect_sizes[(i, j)]

# ------------------------------------------------------------------------
# SINGLE HEATMAP WITH SIGNIFICANCE ASTERISKS
# ------------------------------------------------------------------------
# 1. Create a boolean mask for significance (BH-corrected p < 0.05).
sig_mask = (pval_mat <= 0.05)

# 2. Build an annotation matrix with '*' where significant, and '' otherwise.
#    If you want to ALSO include numeric effect size in each cell, you can
#    do something like: ann_text = f"{effect_mat.loc[i,j]:.2f}*"
annotation_array = np.where(sig_mask, "*", "")


#Computing summary stats
summary = compute_summary_matrices(interaction_values, cell_types, group_labels=(SUBTYPE1, SUBTYPE2))


# 3. Plot the effect_mat as a single heatmap, with asterisks in significant cells.
plt.figure(figsize=(10, 8))
ax = sns.heatmap(
    effect_mat,
    cmap='coolwarm',
    center=0,
    vmin=-1,
    vmax=1,
    annot=annotation_array,  # <-- place asterisks here
    fmt="",                  # <-- no numeric formatting, just the string
    cbar_kws={"label": "Rank-Biserial Correlation"}
)
ax.set_title("Effect Size (Rank-Biserial) with Significant Differences Marked by *")

# Customize x-tick and y-tick labels
ax.set_xticklabels(
    ax.get_xticklabels(),
    fontsize=12,
    #fontweight='bold',
    color='black',
    rotation=45,  # Optional: rotate for better visibility
    ha='right'    # Align labels to the right
)
ax.set_yticklabels(
    ax.get_yticklabels(),
    fontsize=12,
    #fontweight='bold',
    color='black'
)

PLOTS_EXCEL = f'{OUTPUT_DIR}/EXCEL_FILES/'
PLOTS_PNG = f'{OUTPUT_DIR}/PNG/'
PLOTS_SVG = f'{OUTPUT_DIR}/SVG/'

mkdir(PLOTS_EXCEL)
mkdir(PLOTS_PNG)
mkdir(PLOTS_SVG)

plt.tight_layout()

plt.savefig(f'{PLOTS_SVG}/Delta_{SUBTYPE1}_{SUBTYPE2}_ITS_{ITS}.svg', dpi=600, format = 'svg')
plt.savefig(f'{PLOTS_PNG}/Delta_{SUBTYPE1}_{SUBTYPE2}_ITS_{ITS}.png', dpi=600, format='png')

# Define the output Excel file name
output_excel_file = os.path.join(PLOTS_EXCEL, f'Summary_Stats_{SUBTYPE1}_{SUBTYPE2}_ITS_{ITS}.xlsx')

with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
    # -----------------------
    # Sheet 1: Effect Matrix
    # -----------------------
    effect_mat.to_excel(writer, sheet_name='Effect Matrix')

    # -----------------------
    # Sheet 2: P-values From Manwintey U non-parametric test 
    # -----------------------

    pval_mat.to_excel(writer, sheet_name='Multihypotheis Corrected Pvals (Mann-Whitney U test)')
    
    # -----------------------
    # Sheet 3: Annotations
    # -----------------------
    # Convert annotation_array into a DataFrame using effect_mat's index/columns
    annot_df = pd.DataFrame(annotation_array, index=effect_mat.index, columns=effect_mat.columns)
    annot_df = annot_df.astype(str)  # Ensure values are strings
    annot_df.to_excel(writer, sheet_name='Annotations (Significant Only)')
    
    # Force text formatting on annotation cells so * signs are preserved
    ws_annot = writer.sheets['Annotations (Significant Only)']
    # Assuming row 1 and column 1 contain headers/index labels, format cells starting at row 2, col 2:
    for row in ws_annot.iter_rows(min_row=2, max_row=ws_annot.max_row, min_col=2, max_col=ws_annot.max_column):
        for cell in row:
            cell.number_format = '@'
    
    # -----------------------
    # Sheet 4: Mean
    # -----------------------
    mean_sheet_name = 'Mean'
    # Write Subtype1 Mean starting at row 2 (leaving room for a header in row 1)
    summary[SUBTYPE1]['mean'].to_excel(writer, sheet_name=mean_sheet_name, startrow=2, startcol=0)
    # Determine number of rows written (pandas writes a header row, so add 1)
    n_rows_subtype1 = summary[SUBTYPE1]['mean'].shape[0] + 1
    # Write Subtype2 Mean starting a couple of rows below Subtype1 (leave one blank row)
    startrow_subtype2 = n_rows_subtype1 + 2
    summary[SUBTYPE2]['mean'].to_excel(writer, sheet_name=mean_sheet_name, startrow=startrow_subtype2, startcol=0)
    
    # Access the Mean worksheet to add merged header cells
    ws_mean = writer.sheets[mean_sheet_name]
    # Total columns including the index column
    n_cols = summary[SUBTYPE1]['mean'].shape[1] + 1
    # Merge cells in row 1 for the Subtype1 header
    ws_mean.merge_cells(start_row=1, start_column=1, end_row=1, end_column=n_cols)
    ws_mean.cell(row=1, column=1, value=f'{SUBTYPE1} Mean')
    # Merge cells in the row where Subtype2 header should go
    ws_mean.merge_cells(start_row=startrow_subtype2, start_column=1, end_row=startrow_subtype2, end_column=n_cols)
    ws_mean.cell(row=startrow_subtype2, column=1, value=f'{SUBTYPE2} Mean')
    
    # -----------------------
    # Sheet 5: Median
    # -----------------------
    median_sheet_name = 'Median'
    summary[SUBTYPE1]['median'].to_excel(writer, sheet_name=median_sheet_name, startrow=2, startcol=0)
    n_rows_subtype1_med = summary[SUBTYPE1]['median'].shape[0] + 1
    startrow_subtype2_med = n_rows_subtype1_med + 2
    summary[SUBTYPE2]['median'].to_excel(writer, sheet_name=median_sheet_name, startrow=startrow_subtype2_med, startcol=0)
    
    ws_median = writer.sheets[median_sheet_name]
    n_cols_med = summary[SUBTYPE1]['median'].shape[1] + 1
    ws_median.merge_cells(start_row=1, start_column=1, end_row=1, end_column=n_cols_med)
    ws_median.cell(row=1, column=1, value=f'{SUBTYPE1} Median')
    ws_median.merge_cells(start_row=startrow_subtype2_med, start_column=1, end_row=startrow_subtype2_med, end_column=n_cols_med)
    ws_median.cell(row=startrow_subtype2_med, column=1, value=f'{SUBTYPE2} Median')

print(f"Excel file with combined results saved to: {output_excel_file}")