

#### load libraries & utility function 
library(Seurat)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
raw_object_path <- snakemake@input[["raw_object"]]

# outputs
filtered_object_path <- snakemake@output[["filtered_object"]]

# parameters
result_dir <- dirname(filtered_object_path)
filter_expression <- snakemake@params[["filter_expression"]] # "hto_demux != 'Doublet' & hto_demux != 'Negative' & pass_QC == 'True' & guide_call!='Negative' & guide_call!='Multiplet'"
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
                   result_dir=dirname(file.path(filtered_object_path)),
                   prefix='FILTERED_'
                  )