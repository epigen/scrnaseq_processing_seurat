# scRNA-seq Data Processing Snakemake Workflow powered by Seurat

A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for processing and visualizing (multimodal) scRNA-seq data generated with [10X Genomics Kits](https://www.10xgenomics.com/) powered by the R package [Seurat](https://satijalab.org/seurat/index.html).

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)

# Authors
- [Stephan Reichl](https://github.com/sreichl)


# Software
This project wouldn't be possible without the following software

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| Seurat         | https://doi.org/10.1016/j.cell.2021.04.048        |
| SCTransform    | https://doi.org/10.1186/s13059-019-1874-1         |
| inspectdf      | https://github.com/alastairrushworth/inspectdf/   |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

The outlined analyses were performed using the R package Seurat (ver) [ref] unless stated otherwise.

**Merge.** The preprocessed samples were merged using the function Seurat::merge that concatenates the individual samples and their metadata into one Seurat object.

**Metadata.** Metadata was extended with Seurat::PercentageFeatureSet using [X] and by recombination of existing metadata rules [X].

**Guide RNA assignment.** The guide RNA (gRNA) assignment was performed using the protospacer call information provided by the CRISPR functionality of 10x Genomics Cell Ranger (ver). To avoid combinatorial knockout effects and focus on single-gene knockouts only cells with exactly one called protospacer were annotated accordingly, others were labeled as multiplets. 

**Split.** The merged data set was split into subsets by the metadata column(s) [X].

**Filtering.** The cells were filtered by the [X], which resulted in [X] high-quality cells with confident condition and gRNA assignment.

**Normalization.** Filtered count data was normalized using Seurat::SCTransform [ref] with the method parameter glmGamPoi to increase computational efficiency. Other modalities [X] were normalized with Seurat::NormalizeData using method CLR  (Centered Log-Ratio) and margin 2.

**Cell Cycle Scoring.** Cell-cycle scores were determined using the function Seurat::CellCycleScoring using gene lists for M and G2M phase provided by Seurat::cc.genes (Tirosh et al 2015) or [gene lists].

**Cell Scoring.** Cell-module scores were determined using the function Seurat::AddModuleScore using provided gene lists [gene lists].

**Correction.** Filtered count data was first normalized using Seurat::SCTransform [ref] with the method parameter glmGamPoi to increase computational efficiency and identified confounders [X] as covariates to be regressed out.

**Visualization.** To visualize the metadata after each processing step inspectdf (ver) [ref] was used. For the visualization of expression and multimodal [X] data the Seurat functions RidgePlot for ridge-plots, VlnPlot for violin-plots, DotPlot for dot-plots and DoHeatmap for heatmaps were used. Metadata like module scores, were also visualized using RidgePlot for ridge-plots and VlnPlot for violin-plots.

**The processing and analysis described here was performed using a publicly available Snakemake [ver] (ref) workflow [ref - cite this workflow here].**

# Features
The workflow perfroms the following steps. Outputs are always indicated by the respective prefix.
- Preparation (prefix: RAW)
  - Load (mutlimodal) data per 10X sample from Cell Ranger output
  - Load metadata
  - Extend metadata with Seurat::PercentageFeatureSet
  - Assign guide-RNA and KO target gene, if applicable
  - Add metadata columns based on config
  - Create and save Seurat object per sample
- Merge & Split into subsets (prefix: RAW)
  - Merge all samples into one large object, including metadata, called "merged"
  - Split into subsets according to configuration

The following steps are performed on each data split separately.

- Filtering (prefix: FILTERED)
  - Filter cells by a combination of logical-expressions using the metadata
- Normalization (prefix: NORMALIZED)
  - Normalization of expression data using SCTransform, returning normalized values for all expressed genes
  - Normalization of multimodal data using Centered Log-Ratio (CLR)
  - Dynamic highly variable gene (HVG) determination using a residual variance threshold of 1.3 (default)
- Cell Scoring
  - Cell cycle scoring using Seurat::CellCycleScoring and provided S2- and G2M genes
  - Gene module scoring using Seurat::AddModuleScore and provided gene lists
- Correction (prefix: CORRECTED)
  - Normalization and correction for the list of provieded confounders using SCTransform, returning scaled values for all expressed genes
- Visualization by Ridge-, Violin-, Dot-plots and Heatmaps
  - Gene- and feature-lists for plotting can be provided
  - Metadata
    - after each step using inspectdf for all numerical and categorical data
    - ridge- and violin-plots if configured (eg useful for module scores)
  - All expression (RNA) data is plotted after normalization (slot="data") and correction (slot="scale.data"), respectively
  - Ridge plots for other modalities are only generated after normalization (using slot="data")
  - Violin plots for other modalities are only generated after normalization (using slot="data")
  - Dot plots only of normalized data (slot="data"). Performs averaging and scaling (Min-Maxed based) of the values before plotting.
  - Heatmaps always on scaled data (slot="scale.data"). Top 100 HVG are always plotted.
- Reporting (split into categories)
  - Configuration: input configuration, metadata and gene-list files
  - Software: conda environment export to document the installed and used software including versions and build
  - processing_{project_name}: metadata CSV files and visualizations of all processing steps by data split
  - visualization_{project_name}: all supported visualizations by data split
- Save counts
  - functionality to save all counts should be saved as CSV after each processing steps for of all modalities. Useful for downstream analyses that are incompatible with Seurat.
- Results
  - all results will be saved in the result_path as configured where for each data(sub)set a directory with the following structure is created:
  -  counts (for all the .rds object files and .CSV files)
     -  plots (for all visualizations)
     -  stats (for all metadata derived statistics)

# Usage
Here are some tips for the usage of this workflow:
- when generating the sample_annotation sheet use short sample names (they will be the prefix for each barcode in the merged & split datasets)
- run the workflow for each step of processing (with the stop_after parameter) and investigate the results (eg using the report function)
- start with a low complexity in the configuration
- try to finish the analysis of the "merged" data set and later split the data by using the split_by parameter

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
--- COMING SOON ---
