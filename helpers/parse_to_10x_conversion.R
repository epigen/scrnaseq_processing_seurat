# code adapted from
# https://support.parsebiosciences.com/hc/en-us/articles/360053078092-Seurat-Tutorial-65k-PBMCs
# https://rdrr.io/github/MarioniLab/DropletUtils/src/R/write10xCounts.R

# libraries
library(Seurat)
library(dplyr)
library(Matrix)

# config
data_path <- "/path/to/PARSE/data/DGE_filtered"
result_path <- "/path/to/result/directory"
dir.create(result_path, recursive = TRUE)

#### LOAD DATA

# load PARSE count data
data <- ReadMtx(
        mtx = file.path(data_path, "count_matrix.mtx"),
        cells = file.path(data_path, "cell_metadata.csv"),
        features = file.path(data_path, "all_genes.csv"),
        cell.column = 1,
        feature.column = 2,
        cell.sep = ",",
        feature.sep = ",",
        skip.cell = 1,
        skip.feature = 1,
        mtx.transpose = TRUE,
        unique.features = TRUE,
        strip.suffix = FALSE
    )

# load cell metadata from PARSE
cell_meta <- read.csv(file.path(data_path, "cell_metadata.csv"), row.names = 1)

#### TRANSFORM DATA

# check to see if empty gene names are present, add name if so.
table(rownames(data) == "")
rownames(mat)[rownames(data) == ""] <- "unknown"

# transform metadata into desired format (e.g., split or add columns)
# <ADD YOUR CODE HERE>

# create pre-filtered Seurat object to reduce size
data_object <- CreateSeuratObject(data, min.genes = 100, min.cells = 100, names.field = 0, meta.data = cell_meta, project="project_name")

# check the created Seurat object
print(data_object)

#### SAVE DATA

# save Seurat object as RData object
saveRDS(data_object, file=file.path(result_path, "seurat_object.rds"))

# save metadata
write.csv(data_object@meta.data, file.path(result_path, "metadata.csv"))

# save RNA counts
counts_RNA <- data_object@assays$RNA@counts
writeMM(counts_RNA, file=file.path(result_path, "matrix.mtx"))

# save barcodes (i.e., cells)
write(colnames(counts_RNA), file=file.path(result_path, "barcodes.tsv"))

# save features (i.e., genes)
gene.info <- data.frame(rownames(counts_RNA), rownames(counts_RNA), stringsAsFactors=FALSE)
gene.info$gene.type <- rep("Gene Expression", length.out=nrow(gene.info))
write.table(gene.info, file=file.path(result_path, "features.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")