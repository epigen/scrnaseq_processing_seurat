# scRNA-seq Data Processing Snakemake Workflow powered by Seurat

A Snakemake workflow for processing and visualizing (multimodal) scRNA-seq data generated with [10X Genomics Kits](https://www.10xgenomics.com/) powered by the R package [Seurat](https://satijalab.org/seurat/index.html).

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Installation](#installation)
  * [Configuration](#configuration)
  * [Execution](#execution)
  * [Report](#report)
  * [Results](#results)
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

---COMING SOON---

# Features
- Preparation
- Merge & Split into subsets
- Filtering
- Normalization
- Correction
- Visualization
- Reporting



