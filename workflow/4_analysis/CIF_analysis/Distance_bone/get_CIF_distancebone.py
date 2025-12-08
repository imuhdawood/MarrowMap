import os
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from matplotlib import pyplot as plt

# ---------- USER SETTINGS ----------
XL_PATH = 'path/to/cif_distancebone.csv'
KEY_XL_PATH = 'path/to/xenium_keys.csv'
SAVE_PATH = 'path/to/output'
DIST_TO_MICRON = 20 * 0.221

# Visual settings
POINT_ALPHA = 0.6
POINT_SIZE = 40
PALETTE_NAME = 'tab10'
SKIP_MPN = {'PrePMF'}

# Font scaling base: adjust if you prefer bigger/smaller text globally
BASE_FONTSIZE = 12
FIGSIZE = (10, 8)
# -----------------------------------

os.makedirs(SAVE_PATH, exist_ok=True)

# auto-scale fonts with figure dimensions
scale_factor = (FIGSIZE[0] * FIGSIZE[1]) / (10 * 8)
plt.rcParams.update({
    'font.size': BASE_FONTSIZE * scale_factor,
    'axes.titlesize': BASE_FONTSIZE * 1.8 * scale_factor,
    'axes.labelsize': BASE_FONTSIZE * 1.5 * scale_factor,
    'xtick.labelsize': BASE_FONTSIZE * 1.2 * scale_factor,
    'ytick.labelsize': BASE_FONTSIZE * 1.2 * scale_factor,
    'legend.fontsize': BASE_FONTSIZE * 1.2 * scale_factor,
})

key_xl = pd.read_csv(KEY_XL_PATH)
MPN_grp = np.unique(key_xl['Diagnosis_2'])

mpn_list_for_palette = [m for m in MPN_grp if m not in SKIP_MPN]
n_mpn = len(mpn_list_for_palette)
mpn_palette = sns.color_palette(PALETTE_NAME, n_colors=min(10, n_mpn))
if n_mpn > 10:
    mpn_palette = sns.color_palette(None, n_colors=n_mpn)
mpn_color_map = {mpn_list_for_palette[i]: mpn_palette[i % len(mpn_palette)] for i in range(n_mpn)}

min_dist = (((512/20)**2)*2)**0.5
initial_min_dist = ((((128/2)/20)**2)*2)**0.5

for mpn in MPN_grp:
    if mpn in SKIP_MPN:
        continue

    mpn_color = mpn_color_map.get(mpn, (0.2, 0.2, 0.2))
    mpn_idx = np.where(np.array(key_xl['Diagnosis_2']) == mpn)[0]
    xl_list = [
        f"{np.array(key_xl['Slide_No'])[i]}_R{np.array(key_xl['Xenium_region'])[i].split('Region')[-1]}"
        for i in mpn_idx
    ]

    # ---------- ITS-based scatter ----------
    its_group_dist, its_group_cif = [], []
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for sample_id in xl_list:
        csv = os.path.join(XL_PATH, f"{sample_id}.csv")
        if not os.path.exists(csv):
            continue
        df = pd.read_csv(csv)
        if not {'CIF_Score', 'Distance', 'ITS'}.issubset(df.columns):
            continue

        cif, dist, its = df['CIF_Score'].values, df['Distance'].values, df['ITS'].values
        avg_cif, avg_dist = [], []
        its_unq = np.unique(its)
        dist_cnt = int(max(dist) / min_dist) + 1
        for dd in range(dist_cnt):
            for it in its_unq:
                mask_its = its == it
                cif_its, dist_its = cif[mask_its], dist[mask_its]
                if dd == 0:
                    mask = (dist_its >= initial_min_dist) & (dist_its < min_dist)
                else:
                    mask = (dist_its >= min_dist*dd) & (dist_its < min_dist*(dd+1))
                if np.any(mask):
                    avg_cif.append(np.mean(cif_its[mask]))
                    avg_dist.append(np.mean(dist_its[mask]))
        if avg_cif:
            avg_dist = np.array(avg_dist) * DIST_TO_MICRON
            avg_cif = np.array(avg_cif)
            ax.scatter(avg_dist, avg_cif, s=POINT_SIZE, alpha=POINT_ALPHA, color=mpn_color)
            its_group_dist.append(avg_dist)
            its_group_cif.append(avg_cif)

    if its_group_dist:
        gx = np.concatenate(its_group_dist)
        gy = np.concatenate(its_group_cif)
        if len(gx) >= 2:
            r, p = scipy.stats.pearsonr(gx, gy)
            slope, intercept, *_ = scipy.stats.linregress(gx, gy)
            xline = np.linspace(gx.min(), gx.max(), 400)
            yline = intercept + slope * xline
            ax.plot(xline, yline, color='black', linewidth=2.6, label=f"Pearson={r:.3f}, p={p:.3f}")
        ax.set_title(f"{mpn} (ITS-averaged)")
        ax.set_xlabel("Distance to Bone (µm)")
        ax.set_ylabel("CIF")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(SAVE_PATH, f"{mpn}_ITS_scatter.png"), dpi=300)
    plt.close(fig)

print("Done for all MPN groups.")
