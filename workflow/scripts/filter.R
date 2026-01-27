

#### load libraries & utility function 
library("Seurat")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@input[["utils_path"]])

# inputs
raw_object_path <- snakemake@input[["raw_object"]]

# outputs
filtered_object_path <- snakemake@output[["filtered_object"]]

# parameters
filter_expression <- snakemake@params[["filter_expression"]]
# 'flags' for modalities
ab_flag <- snakemake@params[["ab_flag"]]#'AB'
crispr_flag <- snakemake@params[["crispr_flag"]]#'gRNA'
custom_flag <- snakemake@params[["custom_flag"]]#'HTO'

### load raw data
raw_object <- readRDS(file = file.path(raw_object_path))
metadata <- raw_object[[]]

### filter data
if (filter_expression!=""){
    filter_expr <- parse(text=filter_expression)

    filtered_metadata <- subset(metadata, subset= eval(filter_expr))
    filtered_object <- raw_object[,rownames(filtered_metadata)]
}else{
    filtered_metadata <- metadata
    filtered_object <- raw_object
}


### save data
save_seurat_object(seurat_obj=filtered_object,
                   result_dir=dirname(file.path(filtered_object_path))
                  )