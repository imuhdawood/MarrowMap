from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class STConfig:

    # Microns per pixel for spatial scaling
    MPP: float = 0.2125

    # Adipocyte buffer radius in pixels
    ADIPO_BUFFER_RADIUS: int = 15

    # Max allowed distance cell→bone for filtering
    CELL_BONE_MAX_DISTANCE: int = 20000

    # Distance thresholds defining peri-trabecular region (pixels)
    PERITRABECULAR_DISTANCE_THRESHOLD: tuple = (0, 300)

    # Distance threshold defining endosteal region (<1 pixel from bone)
    ENDOSTEAL_DISTANCE_THRESHOLD: int = 1
    
    # Metadata table with Xenium sample keys
    pth_meta_csv: Path = Path("")

    # Output directory for Zarr or analysis artifacts
    pth_out_zar: Path = Path("")

    # Directory containing annotated Xenium cell annotation files
    pth_annotation: Path = Path("")

    # Mapping file linking features → cell types or annotation categories
    pth_feature_map: Path = Path("")

    # CIF (cell interaction feature) scores / coordinates
    pth_cif_score: Path = Path("")

    # Directory to save transformed **negative** annotation results
    pth_neg_annotation_save: Path = Path("")

    # Directory to save transformed **positive** annotation results
    pth_pos_annotation_save: Path = Path("")

    # H&E image directory (OME-TIFF)
    pth_he_img: Path = Path("")

    # Landmark matrix for registration / alignment
    pth_landmark_matrix: Path = Path("")

    # Slimmed metadata table with cell-level annotations
    pth_cell_annotations: Path = Path("")

    # Bone segmentation results directory
    pth_bone_segmentation: Path = Path("")

    # Adipocyte segmentation results directory
    pth_adipocytes_segmentation: Path = Path("")

    # Consolidated AnnData for downstream analysis
    pth_consol_adata: Path = Path("")

    # Output directory for downstream analysis results
    pth_downstream_out: Path = Path("")

    # Directory containing individual sample sdata objects
    pth_sdata: Path = Path("")

    # Final QC’d cell annotation CSV
    pth_cell_annotations_final: Path = Path("")

    # Cell Charter / consensus cell network AnnData
    pth_cell_charter_cn: Path = Path("")

    # Graph-based AnnData for spatial graph modeling
    pth_graph_adata: Path = Path("")

    # Output dir for graph construction results
    pth_graphs_out_dir: Path = Path("")

    # Spatiotype label file (label only)
    pth_spatiotypes_label: Path = Path("")

    # Spatiotype features only (no labels)
    pth_spatiotypes_feat_only: Path = Path("")

    # Spatiotype feature → label mapping file
    pth_spatiotypes_feat_label: Path = Path("")

    # ITS region polygons
    pth_its_poly: Path = Path("")

    # ITS-level spatiotype score predictions
    pth_its_score_file: Path = Path("")

    # Sample-level aggregated spatiotype predictions
    pth_sample_level_score_file: Path = Path("")

    # Bootstrap sample-level predictions for uncertainty estimates
    pth_bootstrap_score_file: Path = Path("")

    # CIF score file (cell-level)
    pth_cif_score_file: Path = Path("")
