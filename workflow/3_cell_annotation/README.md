# Cell Annotation

This segment of the repository contains four sections. Following is the chronological order.

1) Processing data in Python along with all pre-processed elements: `/workflow/3_cell_annotation/1_make_sdata.py`

2) Converting the processed data to RDS objects for use in Seurat in R: `/workflow/3_cell_annotation/2_data_conversion/`
   (First run `workflow/3_cell_annotation/2_data_conversion/2a_zarr_to_anndata.ipynb`, then run `workflow/3_cell_annotation/2_data_conversion/2b_anndata_to_seurat.R`.)

3) Seurat code: `workflow/3_cell_annotation/3_seurat_pipeline.R`

4) Representative plots - UMAPs and heatmaps: `workflow/3_cell_annotation/4_figures.Rmd`, which is wrapped into an HTML here `workflow/3_cell_annotation/4_figures.html`.

## Section 1: Processing in Python
Code: `/workflow/3_cell_annotation/1_make_sdata.py`

We have the following data:
- Raw Xenium files
- CIF scores by tile size 512 and 1024 pixels
- H/E DAPI images
- Positive and negative annotations based on tissue architectural/cytomorphological preservation
- AI-enhanced bone segmentation
- AI-enhanced megakaryocyte segmentation
- Adipocyte segmentation from MarrowQuant 2.0

This script loads the raw Xenium data into AnnData objects. Then it adds in the resegmented megakaryocytes, and aggregates the transcripts based on the new cell boundaries. The rest of the data elements, including H/E images, are added one by one, and the comprehensive data for each sample is stored as a Zarr data object.

## Section 2: Data Format Conversion
Code Folder: `/workflow/3_cell_annotation/2_data_conversion/`

The goal of this section is to convert the Zarr objects made in the previous section into usable Seurat objects, which can then be loaded into R and analyzed using Seurat.
Since there is no direct functionality to do this (at the time of writing this documentation), we convert Zarr -> AnnData -> Seurat object.

Subsection a: `workflow/3_cell_annotation/2_data_conversion/2a_zarr_to_anndata.ipynb`
This section loads the Zarr object, converts it into an AnnData object, and saves it as an .h5ad file.

Subsection b: `workflow/3_cell_annotation/2_data_conversion/2b_anndata_to_seurat.R`
This section loads the AnnData .h5ad file created above and converts it into a Seurat object. This object is ready for downstream processing and analysis.

## Section 3: Seurat Processing in R
Code: `workflow/3_cell_annotation/3_seurat_pipeline.R`

This loads the Seurat-friendly data, QC's it, normalizes using SCTransform, performs linear dimensionality reduction using PCA, runs Harmony for batch-correction and integration, creates a neighborhood graph, clusters using the Leiden algorithm, identifies gene markers for cell type identification/annotation, and runs SingleR for automated cell type annotation. For the latter, the code provided includes an exemplar using the single-cell Azimuth Human Atlas Bone Marrow reference dataset.

## Section 4: Representative Plots
Code Notebook: `workflow/3_cell_annotation/4_figures.Rmd`, which is wrapped into an HTML here `workflow/3_cell_annotation/4_figures.html`.

This code contains exemplars of UMAPs and heatmaps generated to aid cell type annotation.
