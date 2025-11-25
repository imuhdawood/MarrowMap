# config.R
# Dependencies: install.packages("fs") and install.packages("reticulate") if you use python interop
library(fs)

create_STConfig <- function() {
  cfg <- new.env(parent = emptyenv())

  cfg$pth_its_score_file <- path("C:/Xenium/data/NodeProba/CellMicroEnviron/obj_w_Cell_annotations_RC_22122024_final_with_osteo/MIL/ITS_LEVEL/ENVIRON_TWEAK/csv/itr_pred.csv")
  cfg$pth_spatiotypes_its_abundance <- path("C:/Xenium/data/NodeProba/CellMicroEnviron/obj_w_Cell_annotations_RC_22122024_final/MIL/ITS_LEVEL/ENVIRON/Stats.csv")
  cfg$pth_spatiotypes_feat_label <- path('C:/Xenium/data/NodeProba/CellMicroEnviron/obj_w_Cell_annotations_RC_22122024_final_with_osteo/ENVIRON_R2_C3_WITH_FEATS.csv')
  cfg$pth_cell_annotations_final <- path("C:/Xenium/data/NodeProba/CellMicroEnviron/obj_w_Cell_annotations_RC_22122024_final/obj_w_Cell_annotations_RC_22122024_final.csv")
  cfg$seurat_object_rds <- path("C:/xenium_paper/data/july_qc_merged_rds_filtered_3_remove_artefacts_sct_pca_harmony_umap_clus_ds.rds")
  cfg$annotation_with_stromal <- path("C:/xenium_paper/data/final_annotations_extra_cells_with_stromal.csv")
  # Freeze: prevent accidental modification
  lockEnvironment(cfg, bindings = TRUE)
  cfg
}

STConfig <- create_STConfig()
print.STConfig <- function(cfg = STConfig) {
  cat("STConfig (read-only)\n")
  keys <- ls(cfg)
  for (k in keys) {
    val <- get(k, envir = cfg)
    cat(sprintf(" %s: %s\n", k, substr(as.character(val), 1, 200)))
  }
  invisible(cfg)
}
