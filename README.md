# MarrowMap: Redefining the topology of the human bone marrow using augmented spatial transcriptomic analysis

This repository contains code scripts used in the study: **"Redefining the topology of the human bone marrow using augmented spatial transcriptomic analysis"** The study uses image-based spatial transcriptomics data coupled with AI-driven features (structural componets) to define the topology of the normal bone marrow (BM) and BM in myeloproliferative neoplasms (MPNs). A high-level overview of the cohort and analyses conducted is provided in the study workflow figure below. 

![Workflow Diagram](./images/graphical_abstract.png)


## Repository Contents
- **Code Scripts:** Scripts used for data loading, preprocessing, and analysis. 
- **Analysis Outputs:** Where possible we showed the Figures reported in the Manuscript in R or Python Notebooks.

### Usage

#### Step 1: Read Xenium runs data as spatial data object
````python
/workflow/1_xenium_to_sdata/xen_to_sdata.py
````
#### Step 2: Quality Control (removing spatially or morphologically unpreserved, low-transcript cells)
Run the following notebooks
````python
/workflow/2_quality_control/artefact_filtering.ipynb
/workflow/2_quality_control/spatial_preservation.ipynb
````
#### Step 3: Cell annotation
Run the R scripts in the folder `/workflow/3_cell_annotation`. Details about script order is provided in the help file inside the folder.

#### Step 4: Analyses

##### Identifying cell neighbourhgood
- Construct graph representation by running the following script:
````python
workflow/4_analysis/cell_charter/graph_construction.py
````
- Validating graph connectivity and comparision with other graph construction approaches
 ````python
workflow/4_analysis/cell_charter/graph_construction_approaches_comparisions.ipynb
````
##### Abundance Analysis
Run the following notebooks
 ````python
/workflow/4_analysis/cell_lineages_abundance/abundance_plot_fig1.ipynb
/workflow/4_analysis/cell_lineages_abundance/cell_types_viz_trephine.ipynb
````

##### Spatial Analysis (Point Pattern)
Point pattern analysis is done based on the location of corresponding cell/structure and the distance map of the targets. As an example, we will do celltype-celltype analysis.
Run the following scripts for metadata and distance map
````python
/workflow/4_analysis/point_pattern_analysis/get_metadata/get_metadata.py
/workflow/4_analysis/point_pattern_analysis/distance_map/get_distancemap_celltype.py # change the python code depending on which distance map required
````
For structure-structure analysis, we need to compute boundary of structure instead of centre because the size of the structure is much larger compared to the cell.
Run the following script to get the boundary points of the corresponding structure
````python
/workflow/4_analysis/point_pattern_analysis/structure_boundary/get_structure_boundary_pts.py
````
Run the following R script to get point pattern proximity 
````
/workflow/4_analysis/point_pattern_analysis/ppm_Rcode/ppm_celltype_celltype.R # change the R code depending on which analysis you want to compute
````
To plot proximity heatmap we need to gather the corresponding proximity results from R code, compute permutation for p-value, and then plot.
Run following python scripts in order
````python
/workflow/4_analysis/point_pattern_analysis/plot/gather_data/gather_ppm_celltype_celltype.py # change the python code depending on which ppm data you want to gather
# For structure-strucuter make sure run /workflow/4_analysis/point_pattern_analysis/plot/for_permutation first before computing permutation.
/workflow/4_analysis/point_pattern_analysis/plot/for_permutation/get_permutation_pvalue.py
# Note: after computing permutation, create excel (xlsx) file and input the median (from gather_data codes) with the header 'Proximity Rank' and permutation result (from get_permutation_pvalue.py) with the header 'P-Value'.
/workflow/4_analysis/point_pattern_analysis/plot/proximity_heatmap/get_proximity_heatmap.py
````
Run the following python script to do delta analysis with the corresponding heatmap
````python
/workflow/4_analysis/point_pattern_analysis/plot/for_delta/get_delta.py
````
##### Visualization of niches
Run the following notebooks
 ````python
/workflow/4_analysis/spatial_analysis/visualize_niches.ipynb
/workflow/4_analysis/spatial_analysis/interactions_rank.py
````

##### Spatiotypes Analysis
The spatiotypes analysis was done using band descriptor.
Run the following scripts in order
 ````python
/workflow/4_analysis/spatioytypes/compute_band_descriptors.py
/workflow/4_analysis/spatioytypes/cluster_descriptors.py
/workflow/4_analysis/spatioytypes/spatiotypes_viz_trephine.ipynb
````
Run the R notebook to visualize the spatiotypes spatial cellular abundance
````
/workflow/4_analysis/spatioytypes/spatiotypes_defination.Rmd
````

##### Multiple Instance Learning Model
Run the following scripts in order
 ````python
/xenium_paper/workflow/5_MIL/main.py
/workflow/5_MIL/visualization/MIL_figures.ipynb
````
For visualizing intra- and inter- samples hetrogenity using Circos plot

```
/workflow/5_MIL/visualization/circos_plots.Rmd
```

## Citation

If you use this repository, please cite our preprint:

> [Redefining the topology of the human bone marrow using augmented spatial transcriptomic analysis](https://www.biorxiv.org/content/10.1101/2024.06.23.600257v1)

## License

The source code in this repo is released under the MIT-CC-Non-Commercial license. For external packages (CellCharter, BandDescriptor, ESIML) please refer to check corresponding code repositories.

---

For more information, please refer to the [preprint](https://www.biorxiv.org/content/10.1101/2024.06.23.600257v1).
