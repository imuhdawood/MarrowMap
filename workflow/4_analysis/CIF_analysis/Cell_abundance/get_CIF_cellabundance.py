import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
from matplotlib.ticker import MultipleLocator

anno_numb = 3
df_main = pd.read_csv('path/to/main.csv', index_col=0)

# -------------------- USER CONFIG: update these --------------------
xl_path = 'path/to/CIF_cellabundance.xlsx'   # <-- change this to your data folder
save_path = 'path/to/output'
# ------------------------------------------------------------------

# Column names in CSVs (cell types)
fts_grp = np.array(df_main['obj.anno_{}_w_megs_w_stromal'.format(anno_numb)])
fts_grp = fts_grp.astype(str)
abund_idxs = np.unique(fts_grp)
# remove a few specific entries present in original script
abund_idxs = np.delete(abund_idxs,[np.where(abund_idxs=='nan')][0][0])
abund_idxs = np.delete(abund_idxs,[np.where(abund_idxs=='CD69')][0][0])
abund_idxs = np.delete(abund_idxs,[np.where(abund_idxs=='Adipo-MSC')][0][0])
cell_type = abund_idxs.tolist()

# MPN groups (subfolder names inside xl_path)
mpn_grp = ["Normal", "ET", "PV", "MF"]

# Color mapping for MPNs (used in scatter).
colors_mapping = {
    "Normal": '#008000',
    "ET": '#0000FF',
    "PV": '#FF0000',
    "MF": '#FFA500'
}

# Ensure save folder exists
os.makedirs(save_path, exist_ok=True)

# Build the list of MPN groups to actually use, excluding PrePMF
mpn_grp_used = [g for g in mpn_grp if g != "PrePMF"]
if len(mpn_grp) != len(mpn_grp_used):
    print("Info: 'PrePMF' detected in mpn_grp and will be excluded from plotting for now.")

# -------------------- ADDED: helper functions for CI --------------------

def bootstrap_ci(x, y, func, n_boot=2000, alpha=0.05, random_state=None):
    """
    Compute bootstrap percentile CI for a statistic function func(x_sample, y_sample) -> scalar or array.
    func should accept (x_resampled, y_resampled) and return a numpy array (1D).
    Returns: (est, lower, upper, boot_dist)
    """
    rng = np.random.default_rng(random_state)
    n = len(x)
    if n == 0:
        raise ValueError("No data for bootstrap.")
    est = func(x, y)
    boot_samples = []
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            boot_val = func(x[idx], y[idx])
            boot_samples.append(np.asarray(boot_val))
        except Exception:
            boot_samples.append(np.full_like(np.asarray(est), np.nan, dtype=float))
    boot_arr = np.vstack(boot_samples)  # shape (n_boot, len(est))
    lower = np.nanpercentile(boot_arr, 100 * (alpha / 2.0), axis=0)
    upper = np.nanpercentile(boot_arr, 100 * (1 - alpha / 2.0), axis=0)
    return np.asarray(est), np.asarray(lower), np.asarray(upper), boot_arr

def pearson_r(x, y):
    """Return Pearson r (scalar). If invalid, return np.nan"""
    try:
        if len(x) < 2:
            return np.nan
        r, _ = scipy.stats.pearsonr(np.asarray(x).astype(float), np.asarray(y).astype(float))
        return np.asarray(r)
    except Exception:
        return np.nan

def fit_linear_params(x, y):
    """Return slope, intercept as numpy array [slope, intercept] using np.polyfit."""
    try:
        if len(x) < 2:
            return np.array([np.nan, np.nan])
        slope, intercept = np.polyfit(np.asarray(x).astype(float), np.asarray(y).astype(float), 1)
        return np.array([slope, intercept])
    except Exception:
        return np.array([np.nan, np.nan])

def predict_from_params(params, x_vals):
    """params: [slope, intercept]"""
    slope, intercept = params
    return slope * x_vals + intercept

def fit_and_plot_ci(ax, x, y, color='0.1', n_boot=2000, alpha=0.05):
    """
    Fit linear regression and plot point estimate and bootstrap CI band on ax.
    - CI band: grey, semi-transparent, no legend entry
    - fit line: solid, no legend entry
    Returns dictionary with slope/intercept and Pearson r + bootstrap CI arrays.
    """
    x = np.asarray(x).astype(float)
    y = np.asarray(y).astype(float)
    valid = ~np.isnan(x) & ~np.isnan(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 2:
        return {
            'slope': (np.nan, np.nan, np.nan),
            'intercept': (np.nan, np.nan, np.nan),
            'r': (np.nan, np.nan, np.nan),
            'pred_grid': None
        }

    # Prepare x grid for plotting CI band
    x_grid = np.linspace(np.min(x), np.max(x), 200)

    # bootstrap for predicted values at each x_grid point
    def pred_func(xs, ys):
        p = fit_linear_params(xs, ys)
        return predict_from_params(p, x_grid)  # returns array shape (len(x_grid),)

    est_pred, lower_pred, upper_pred, pred_boot = bootstrap_ci(x, y, pred_func, n_boot=n_boot, alpha=alpha)

    # For slope/intercept CIs
    def params_func(xs, ys):
        return fit_linear_params(xs, ys)  # returns array([slope, intercept])
    params_est_full, params_low, params_up, params_boot = bootstrap_ci(x, y, params_func, n_boot=n_boot, alpha=alpha)

    # For Pearson r CI (bootstrap)
    _, r_low, r_up, r_boot = bootstrap_ci(x, y, lambda a,b: pearson_r(a,b), n_boot=n_boot, alpha=alpha)

    # Plot the CI band and fit line — grey shading, no legend entries, solid line
    ax.fill_between(
        x_grid, lower_pred, upper_pred,
        color='grey', alpha=0.25, edgecolor=None  # <-- grey shading, no label
    )
    ax.plot(x_grid, est_pred, linestyle='-', linewidth=1.5, color=color)  # solid line, no label

    return {
        'slope': (params_est_full[0], params_low[0], params_up[0]),
        'intercept': (params_est_full[1], params_low[1], params_up[1]),
        'r': (np.asarray(pearson_r(x, y)), r_low, r_up),
        'pred_grid': (x_grid, est_pred, lower_pred, upper_pred),
        'boots': {
            'params_boot': params_boot,
            'pred_boot': pred_boot,
            'r_boot': r_boot
        }
    }

# -------------------- END helper functions --------------------

# Prepare the overall grid (fig_all)
n_cells = len(cell_type)
cols_all = 4
rows_all = math.ceil(n_cells / cols_all)
fig_all, axes_all = plt.subplots(rows_all, cols_all, figsize=(cols_all * 6, rows_all * 5), squeeze=False)
plt.subplots_adjust(hspace=0.4, wspace=0.35)

cnt_row = 0
cnt_col = 0

for c in range(len(cell_type)):
    ct = cell_type[c]

    # Collect data rows (abundance, CIF, MPN) for this cell type across used MPN groups
    rows = []
    for mpn in mpn_grp_used:
        sub_path = os.path.join(xl_path, mpn)
        if not os.path.isdir(sub_path):
            # folder missing: skip with informative message
            print(f"Warning: folder not found: {sub_path} (skipping {mpn})")
            continue

        csv_files = [f for f in os.listdir(sub_path) if f.lower().endswith(".csv")]
        for fname in csv_files:
            fullp = os.path.join(sub_path, fname)
            try:
                xl = pd.read_csv(fullp)
            except Exception as e:
                print(f"Warning: could not read {fullp}: {e}")
                continue

            if ct not in xl.columns:
                continue

            for ridx in range(len(xl)):
                try:
                    abundance_val = xl.iloc[ridx, xl.columns.get_loc(ct)]
                    cif_val = xl.iloc[ridx, -1]  # assume CIF is last column
                    if pd.isna(abundance_val) or pd.isna(cif_val):
                        continue
                    rows.append((float(abundance_val), float(cif_val), mpn))
                except Exception:
                    continue

    # Build DataFrame for this cell type
    if len(rows) == 0:
        df = pd.DataFrame(columns=[ct, "CIF", "MPN"])
    else:
        df = pd.DataFrame(rows, columns=[ct, "CIF", "MPN"])

    # determine per-cell figure size based on number of samples (so not every cell plot is same size)
    total_samples = len(df)
    # heuristic scaling: base size plus small increase with samples, bounded to reasonable sizes
    fig_w = min(22, max(8, 6 + total_samples * 0.08))
    fig_h = min(18, max(6, 5 + total_samples * 0.04))
    fig_mpn, axes_mpn = plt.subplots(2, 2, figsize=(fig_w, fig_h))
    axes_flat = axes_mpn.flatten()

    for j, mpn in enumerate(mpn_grp_used):
        if j >= 4:
            print(f"Note: more than 4 MPN groups in mpn_grp_used; only first 4 are shown in 2x2. (mpn {mpn} at index {j} is ignored in 2x2 layout)")
            break
        ax = axes_flat[j]
        df_m = df[df["MPN"] == mpn] if not df.empty else pd.DataFrame(columns=[ct, "CIF", "MPN"])

        if df_m.shape[0] == 0:
            ax.text(0.5, 0.5, f'No samples for {mpn}', ha='center', va='center', fontsize=25)
            ax.set_title(f'{mpn}', fontsize=22, fontweight='bold')  # title bold, no n
            ax.set_xlabel(ct, fontsize=22)
            ax.set_ylabel("CIF Score", fontsize=22)
            continue

        # draw points (label only for samples — CI & fit have no labels so they won't appear in legend)
        color = colors_mapping.get(mpn, "k")
        ax.scatter(df_m[ct], df_m["CIF"], s=200, c=color, alpha=0.9, label='samples')

        # seaborn regression (just the underlying fit visualization if desired) without its own CI
        try:
            sns.regplot(x=ct, y="CIF", data=df_m, scatter=False, ax=ax, truncate=True, ci=None, color=".1")
        except Exception as e:
            print(f"Warning: regplot failed for {ct}/{mpn}: {e}")

        # Fit and plot bootstrap CI band + fit line (our added function)
        try:
            res = fit_and_plot_ci(ax, df_m[ct].values, df_m["CIF"].values, color='k', n_boot=1000, alpha=0.05)
        except Exception as e:
            print(f"Warning: CI plotting failed for {ct}/{mpn}: {e}")
            res = None

        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.set_xlabel('Abundance', fontsize=22, fontweight='bold')
        ax.set_ylabel("CIF Score", fontsize=22, fontweight='bold')

        # Pearson correlation (original + display in title)
        r_val = np.nan
        p_val = np.nan
        if df_m.shape[0] >= 2:
            try:
                r_val, p_val = scipy.stats.pearsonr(df_m[ct].astype(float), df_m["CIF"].astype(float))
            except Exception:
                r_val, p_val = np.nan, np.nan

        # Title includes MPN name and Pearson r and p-value with three decimals, bold; no sample counts
        if not np.isnan(r_val):
            title_r = f"r={r_val:.3f} (p={p_val:.10f})"
        else:
            title_r = "r=NA"
        ax.set_title(f'{mpn}: {title_r}', fontsize=22, fontweight='bold')


    # If fewer than 4 used MPN groups, blank remaining axes
    for k in range(len(mpn_grp_used), 4):
        axes_flat[k].axis("off")

    # Suptitle for the per-cell figure is only the cell name (ct), bold
    fig_mpn.suptitle(f'{ct}', fontsize=50, fontweight='bold')
    fig_mpn.tight_layout(rect=[0, 0.03, 1, 0.95])

    safe_name = ct.replace("/", "_")
    outfn_mpn = os.path.join(save_path, f"{safe_name}.png")
    fig_mpn.savefig(outfn_mpn, dpi=300, bbox_inches="tight")
    plt.close(fig_mpn)


    # advance grid index
    if cnt_col == cols_all - 1:
        cnt_col = 0
        cnt_row += 1
    else:
        cnt_col += 1

# hide unused axes in overall grid
for r in range(rows_all):
    for c in range(cols_all):
        idx = r * cols_all + c
        if idx >= n_cells:
            axes_all[r][c].axis("off")


print("Done. Saved per-cell 2x2 figures (excluding PrePMF) and the overall grid to:", save_path)
