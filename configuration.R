# config.R
# Dependencies: install.packages("fs") and install.packages("reticulate") if you use python interop
library(fs)

create_STConfig <- function() {
  cfg <- new.env(parent = emptyenv())

  # Path to ITS-level spatiotype score predictions (per-cell or per-region predictions)
  cfg$pth_its_score_file <- ""

  # Path to ITS-level spatiotype abundance summary statistics
  cfg$pth_spatiotypes_its_abundance <- ""

  # Path to CSV mapping spatiotype features to labels (feature/label reference file)
  cfg$pth_spatiotypes_feat_label <- ""

  # Path to final cell annotation table used throughout the analysis
  cfg$pth_cell_annotations_final <- ""

  # Path to merged and QC-filtered Seurat object (.rds)
  cfg$seurat_object_rds <- ""

  # Path to annotation file including stromal cell annotations
  cfg$annotation_with_stromal <- ""

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
