
#### load libraries & utility function 
library("Seurat")
library("data.table")
library("dplyr")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R")

# inputs
filtered_object_path <- snakemake@input[["filtered_object"]]

# outputs
pseudobulk_counts_path <- snakemake@output[["pseudobulk_counts"]]

# parameters
ab_flag <- snakemake@config[["modality_flags"]][['Antibody_Capture']]
crispr_flag <- snakemake@config[["modality_flags"]][['CRISPR_Guide_Capture']]
custom_flag <- snakemake@config[["modality_flags"]][['Custom']]

pseudobulk_by <- snakemake@config[["pseudobulk"]][["by"]]
pseudobulk_method <- match.fun(snakemake@config[["pseudobulk"]][["method"]]) # can be "sum", "mean", or "median"

### load filtered data
seurat_object <- readRDS(file = file.path(filtered_object_path))

# get metadata
metadata <- as.data.frame(seurat_object[[]])
metadata$ID <- rownames(metadata)

# pseudobulk per modality (RNA, ...)
for (modality in c("RNA", ab_flag, crispr_flag, custom_flag)){
    if (modality==""){
        next
    }

    # extract data
    tmp_data <- as.data.frame(t(as.data.frame(GetAssayData(object = seurat_object, slot = "counts", assay = modality))))
    features <- colnames(tmp_data)
    tmp_data$ID <- rownames(tmp_data)

    # join with metadata
    tmp_data <- inner_join(tmp_data, metadata, by = "ID")

    # pseudobulk by method
    tmp_pseudobulk <- tmp_data %>%
      group_by(across(all_of(pseudobulk_by))) %>%
      summarise(across(all_of(features), pseudobulk_method, .names = "{.col}"), .groups = "drop")

    # reformat df and create metadata
    tmp_pseudobulk <- as.data.frame(tmp_pseudobulk)
    rownames(tmp_pseudobulk) <- apply(tmp_pseudobulk[, pseudobulk_by], 1, function(x) paste(x, collapse = "_"))
    tmp_pseudobulk[,pseudobulk_by] <- NULL
    tmp_pseudobulk <- as.data.frame(t(tmp_pseudobulk))

    # save pseudobulked data frame
    fwrite(as.data.frame(tmp_pseudobulk), file=file.path(dirname(pseudobulk_counts_path),paste0(modality,".csv")), row.names=TRUE)
}


# generate aggregated metadata sheet
metadata_aggregated <- metadata %>%
  group_by(across(all_of(pseudobulk_by))) %>%
  summarise(cell_count = n(), .groups = "drop")
# format metadata
metadata_aggregated <- as.data.frame(metadata_aggregated)
rownames(metadata_aggregated) <- apply(metadata_aggregated[, pseudobulk_by], 1, function(x) paste(x, collapse = "_"))
# save metadata
fwrite(as.data.frame(metadata_aggregated), file=file.path(dirname(pseudobulk_counts_path),"metadata.csv"), row.names=TRUE)