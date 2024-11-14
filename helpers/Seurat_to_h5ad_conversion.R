# adapted from: https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html

# libraries
library("Seurat")
library("SeuratDisk")

# configs
object_path <- file.path("path/to/Seurat/object.rds")
results_path <- file.path("result/path")
dir.create(results_path, recursive = TRUE)

# load seurat object
data <- readRDS(file = file.path(object_path))

# save as h5Seurat
SaveH5Seurat(data, filename = file.path(results_path,"object.h5Seurat"), overwrite = TRUE)

# convert to anndata for scanpy (https://mojaveazure.github.io/seurat-disk/reference/Convert.html)
Convert(file.path(results_path,"object.h5Seurat"), assay = "SCT", dest = "h5ad", overwrite = TRUE)
