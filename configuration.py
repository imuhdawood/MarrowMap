from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class STConfig:
    pth_meta_csv: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/xenium_keys.csv")
    pth_out_zar: Path = Path("/well/rittscher/users/qwi813/xenium_paper/outputs")
    pth_annotation: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/annotations/Xenium_annotations_DAN/")
    pth_cif_score: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/CIF/Coordinates_Scores/")
    pth_neg_annotation_save: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/transformed_negative_annotations/")
    pth_he_img: Path = Path("/well/rittscher/users/qdv200/MPN/xenium/data/morphology_files/ome_tif/")
    pth_landmark_matrix: Path = Path('/well/rittscher/users/qdv200/MPN/xenium/data/landmarks_matrix/')