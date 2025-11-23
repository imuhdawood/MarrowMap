cell_color_dict = { 
'Erythroid': '#EE0000', 
'SMC': '#CD5C5C', 
'GMP': '#9900b3', 
'HSPC': '#d58ad5', 
'Monocyte': '#2384c6', 
'Macrophage': '#5F9EA0', 
'Myeloid': '#b3ddf2', 
'Granulocyte/mast': '#7bbce8', 
'Stromal': '#f3a314', 
'Endothelial': '#ffcc00', 
'Adipo-MSC': '#ffc0cb', 
'T_cell': '#33b233', 
'B_cell': '#b0f2b0', 
'CD69': '#1db388', 
'Plasma_cell': '#7edc89', 
'Osteo-MSC': '#FF7256', 
'Megakaryocyte': '#d633d6', 
'DC': '#00eeee', 
'adipocyte': '#ffdab3',
'Unknown' : "#b6babc" 
} 


cell_color_dict_coarse = {
    "Stromal": "pink",
    "Lymphocyte": "green",
    "HSPC": "brown",
    "Endothelial": "blue",
    "Immature_myeloid": "#00778B",
    "Myeloid": "#E6A800",
    "Erythroid": "red",
    "MNP": "skyblue",
    "Osteo-MSC": "#6A5ACD",
    "Megakaryocyte": "purple",
    "Granulocyte/mast": "gray"
}


cn_color_palette = { 
    "0": "blue", 
    "1": "orange", 
    "2": "green", 
    "3": "red", 
    "4": "purple", 
    "5": "brown", 
    "6": "pink", 
    "7": "gray", 
    "8": "olive", 
    "9": "cyan", 
    "Unknown":"#b6babc" 
} 


# spatiotypes_color_map = {
#   "Stromal-0": "#fed7de",
#   "Stromal-1": "#ffc0cb",
#   "Stromal-2": "#ff3659",
#   "Lymphocyte-0": "#16f916",
#   "Lymphocyte-1": "#008000",
#   "Lymphocyte-2": "#003300",
#   "HSPC-0": "#d56a6a",
#   "HSPC-1": "#a52a2a",
#   "HSPC-2": "#4d1212",
#   "Endothelial-0": "#6060fb",
#   "Endothelial-1": "#0000ff",
#   "Endothelial-2": "#000085",
#   "Immature_myeloid-0": "#21daf9",
#   "Immature_myeloid-1": "#00778b",
#   "Immature_myeloid-2": "#002c33",
#   "Myeloid-0": "#fbd05b",
#   "Myeloid-1": "#e6a800",
#   "Myeloid-2": "#765600",
#   "Erythroid-0": "#fb6060",
#   "Erythroid-1": "#ff0000",
#   "Erythroid-2": "#850000",
#   "MNP-0": "#c5e6f4",
#   "MNP-1": "#87ceeb",
#   "MNP-2": "#1d99cb",
#   "Osteo-MSC-0": "#aba2e1",
#   "Osteo-MSC-1": "#6a5acd",
#   "Osteo-MSC-2": "#332687",
#   "Megakaryocyte-0": "#f916f9",
#   "Megakaryocyte-1": "#800080",
#   "Megakaryocyte-2": "#330033",
#   "Granulocyte/mast-0": "#aeaeae",
#   "Granulocyte/mast-1": "#808080",
#   "Granulocyte/mast-2": "#434343"
# }

spatiotypes_color_map = {
    # --- STROMAL (Rose/Mauve/Bean) ---
    "Stromal-0": "#B05B6F",   # Antique Rose
    "Stromal-1": "#885578",   # Muted Mauve
    "Stromal-2": "#452C2C",   # Dark Bean

    # --- LYMPHOCYTE (Orange/Bronze) ---
    "Lymphocyte-0": "#FF913F", # Carrot Orange
    "Lymphocyte-1": "#D16100", # Burnt Orange (Override Requested)
    "Lymphocyte-2": "#7A4900", # Bronze

    # --- HSPC (Pinks/Maroons) ---
    "HSPC-0": "#FF2F80",      # Hot Pink
    "HSPC-1": "#CC0744",      # Crimson
    "HSPC-2": "#5A0007",      # Deep Maroon

    # --- ENDOTHELIAL (Blues) ---
    "Endothelial-0": "#0AA6D8", # Cerulean Blue
    "Endothelial-1": "#3B5DFF", # Royal Blue
    "Endothelial-2": "#000035", # Midnight Blue

    # --- IMMATURE MYELOID (Teal/Dark Green) ---
    "Immature_myeloid-0": "#00C2A0", # Turquoise
    "Immature_myeloid-1": "#00846F", # Teal
    "Immature_myeloid-2": "#001E09", # Black-Green

    # --- MYELOID (Gold/Dark Browns) ---
    "Myeloid-0": "#FFB500",   # Goldenrod
    "Myeloid-1": "#372101",   # Dark Chocolate
    "Myeloid-2": "#1E0200",   # Black-Red

    # --- ERYTHROID (Reds) ---
    "Erythroid-0": "#FF4A46", # Salmon Red
    "Erythroid-1": "#E83000", # Strong Red-Orange
    "Erythroid-2": "#BA0900", # Brick Red

    # --- MNP (Ocean/Navy/Grey) ---
    "MNP-0": "#006FA6",       # Ocean Blue
    "MNP-1": "#012C58",       # Deep Navy
    "MNP-2": "#61615A",       # Dark Grey

    # --- OSTEO-MSC (Purples) ---
    "Osteo-MSC-0": "#A079BF", # Deep Lavender
    "Osteo-MSC-1": "#72418F", # Purple
    "Osteo-MSC-2": "#320033", # Dark Plum

    # --- MEGAKARYOCYTE (Greens) ---
    "Megakaryocyte-0": "#549E79", # Sage Green (Override: Lightest visible green)
    "Megakaryocyte-1": "#008941", # Standard Green
    "Megakaryocyte-2": "#004D43", # Deep Forest Green

    # --- GRANULOCYTE/MAST (Violet/Olive/Bordeaux) ---
    "Granulocyte/mast-0": "#7900D7", # Vivid Violet
    "Granulocyte/mast-1": "#456648", # Olive Drab
    "Granulocyte/mast-2": "#6B002C", # Bordeaux
}