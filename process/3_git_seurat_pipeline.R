options(future.globals.maxSize = 8000 * 1024^3) #increase font size

#load required libraries
library(Seurat)
library(ggplot2)
library(harmony)
library(reticulate)

#read the merged object into environment
obj = readRDS("<insert path to merged object>")

#filter out cells that are empty (i.e. contain no transcripts.)
obj <- subset(obj, subset = nCount_RNA > 0)

# Run SCTransform
obj <- SCTransform(obj, assay = "RNA")

# Run PCA
obj <- RunPCA(obj, npcs = 50, features = rownames(obj))

# Run Harmony
obj <- RunHarmony(obj,
 				group.by.vars = c("run"),
 				reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

#finding neighbors
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)

#finding clusters
obj <- FindClusters(obj, resolution = 0.5)

#install leiden algorithm if not done already
#reticulate::py_install("leidenalg", forge = TRUE)

#finding neighbors with fewer dimensions
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, k.param = 20, prune.SNN = 0.2)

# Select a range of resolutions
resolution.range <- seq(from = 0.2, to = 1, by = 0.2)
#finding clusters again
obj <- FindClusters(obj, resolution = resolution.range, group.singletons = TRUE)

# identify cell types
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "wilcox", adj.p-val < 0.05 )

# run singleR using azimuth ref dataset

# loading libraries
library(rlang)
library(SingleR)
library(Seurat)
library(dplyr)
library(tidyverse)

azimuth_bm_ref = readRDS("<insert path containing ref RDS object; can be obtained by emailing azimuth>") #reading the reference object containing the azimuth bone marrow reference data

intersect_genes = intersect(rownames(azimuth_bm_ref), rownames(obj)) #finding common genes between our panel and the reference

azimuth_bm_ref <- subset(azimuth_bm_ref, features = intersect_genes) #subsetting ref data to only contain the common genes (since we don't need the others)

azimuth_bm_ref_counts = LayerData(azimuth_bm_ref, assay = "RNA", layer = 'data') #extracting raw counts from ref

obj_counts <- LayerData(obj, assay = "Xenium", layer = 'counts') #extracting raw counts for downsampled data

#running the annotation algorithm 
azimuth_anns <- SingleR(test = obj_counts, 
                  ref = azimuth_bm_ref_counts, 
                  labels = azimuth_bm_ref$celltype.l2,
                  de.method = 'wilcox') #crashed twice, submitting via a SLURM job

#cleaning environment
rm(intersect_genes)
rm(azimuth_bm_ref)
rm(azimuth_bm_ref_counts)
rm(obj_counts)
