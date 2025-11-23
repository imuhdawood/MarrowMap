from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class STConfig:
    pth_meta_csv: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/xenium_keys.csv")
    pth_out_zar: Path = Path("/well/rittscher/users/qwi813/xenium_paper/outputs")
    pth_annotation: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/annotations/Xenium_annotations_DAN/")
    pth_feature_map: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/cell_annotation/revised_feature_map.csv")
    pth_cif_score: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/CIF/Coordinates_Scores/")
    pth_neg_annotation_save: Path = Path("/well/rittscher/users/qwi813/xenium_paper/outputs/transformed_negative_annotations/")
    pth_pos_annotation_save: Path = Path("/well/rittscher/users/qwi813/xenium_paper/outputs/transformed_positive_annotations/")
    pth_he_img: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/morphology_files/ome_tif/")
    pth_landmark_matrix: Path = Path('/well/rittscher/users/qdv200/MPN/xenium/data/landmarks_matrix/')
    pth_cell_annotations: Path = Path('/well/rittscher/users/qdv200/MPN/xenium/data/object_metadata/slimmed_all_metadata_03.csv')
    pth_bone_segmentation: Path = Path('/well/rittscher/users/qdv200/MPN/xenium/data/bone_segmentation/')
    pth_adipocytes_segmentation: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/adipocytes/")
    MPP: float = 0.2125  # microns per pixel
    ADIPO_BUFFER_RADIUS: int = 15  # pixels
    CELL_BONE_MAX_DISTANCE: int = 20000  # pixels
    PERITRABECULAR_DISTANCE_THRESHOLD: tuple = (0, 300)  # pixels
    ENDOSTEAL_DISTANCE_THRESHOLD: int = 1 # less than 1 pixel distance to bone


    # Data paths for downstream analysis
    pth_consol_adata: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/july_sdata_with_adipocytes/all_merged_sdata_with_adipocytes.h5ad")
    pth_downstream_out: Path = Path("/well/rittscher/users/qwi813/xenium_paper/outputs/downstream_analysis/")
    pth_sdata: Path = Path('/well/rittscher/users/qdv200/MPN/xenium/data/july_sdata/')
    pth_cell_annotations_final: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/july_merged_rds/obj_w_Cell_annotations_RC_22122024_final.csv")
    pth_cell_charter_cn: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/paper_code/paper_annotations/data/anndata_cc_3_n10.h5ad")
    pth_graph_adata: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/paper_code/3_analysis/4_graph_creation/data/10693_R1_sp.h5ad")
    pth_graphs_out_dir: Path = Path("/well/rittscher/users/qwi813/xenium_paper/outputs")

    pth_spatiotypes_label: Path = Path("/well/rittscher/users/qwi813/dev/XeniumPipe/xen_paper_plots/data/ENVIRON_R2_C3.csv")