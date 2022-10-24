
#### load libraries & utility function 
library(Seurat)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
merged_object_path <- snakemake@input[["merged_object"]]

# outputs
split_object_path <- snakemake@output[["split_object"]]

# parameters
result_dir <- snakemake@params[["result_dir"]]
split_by <- snakemake@params[["split"]] #"condition_8h_cytokines"
# 'flags' for modalities
ab_flag <- snakemake@params[["ab_flag"]]#'AB'
crispr_flag <- snakemake@params[["crispr_flag"]]#'gRNA'
custom_flag <- snakemake@params[["custom_flag"]]#'HTO'

### load merged data
merged_object <- readRDS(file = file.path(merged_object_path))
metadata <- merged_object[[]]


### split data
split_list <- regmatches(split_by, regexpr("_", split_by), invert = TRUE)
split <- split_list[[1]][1]
cat <- split_list[[1]][2]

split_expr <- parse(text=paste0(split,'==',deparse(cat)))

tmp_metadata <- subset(metadata, subset= eval(split_expr))
tmp_object <- merged_object[,rownames(tmp_metadata)]

### save data
save_seurat_object(seurat_obj=tmp_object,
                   result_dir=dirname(split_object_path),
                   prefix='RAW_'
                  )
