#install anndata if not done already
packageurl <- "https://cran.r-project.org/src/contrib/Archive/anndata/anndata_0.7.5.tar.gz" #10.4
install.packages(packageurl, repos=NULL, type="source")

#loading libraries
library(Seurat)
library(sceasy)
library(reticulate)
library(anndata)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

#finding anndata files
base_path = "<insert path to directory that contains all anndata files; this will be used as working directory>"
files = list.files(base_path) #listing files in abovementioned directory
files = files[grep(".h5ad", files)] #grepping files with the anndata file extension

for (i in files){
  print(i) #logging progress
  file_path = paste0(base_path, i)
  data <- read_h5ad(file_path) #reading anndata into environment 
  expr_matrix <- t(as.matrix(data$X)) #extracting matrix
  rownames(expr_matrix) = data$var_names
  data <- CreateSeuratObject(counts = expr_matrix, meta.data = data$obs) #creating seurat object
  saveRDS(data, sub("h5ad", "rds", file_path)) #saving seurat object
}

#make cell_id the rownames
base_path = "<insert path to working directory; will normally be same as above>"
files = list.files(base_path)
files = files[grep("rds", files)] #subset to just the rds files

for (i in files){
  print(i) #logging progress
  file_path = paste0(base_path, i)
  data <- readRDS(file_path)
  rownames(data@meta.data) <- data$cell_id
  saveRDS(data, sub("rds", "rds", file_path))
}