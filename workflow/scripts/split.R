
#### load libraries & utility function 
library("Seurat")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# inputs
merged_object_path <- snakemake@input[["merged_object"]]

# outputs
split_object_path <- snakemake@output[["split_object"]]

# parameters
split_by <- snakemake@wildcards[["split"]]
# 'flags' for modalities
ab_flag <- snakemake@config[["modality_flags"]][['Antibody_Capture']]
crispr_flag <- snakemake@config[["modality_flags"]][['CRISPR_Guide_Capture']]
custom_flag <- snakemake@config[["modality_flags"]][['Custom']]

### load merged data
merged_object <- readRDS(file = file.path(merged_object_path))
metadata <- merged_object[[]]

if(split_by=="merged"){
    tmp_object <- merged_object
}else{
    ### split data
    split_list <- unlist(regmatches(split_by, regexpr("__", split_by), invert = TRUE))
    split <- split_list[1]
    cat <- split_list[2]
    
    split_expr <- parse(text=paste0(split,'==',deparse(cat)))
    
    tmp_metadata <- subset(metadata, subset= eval(split_expr))
    tmp_object <- merged_object[,rownames(tmp_metadata)]
}

### save data
save_seurat_object(seurat_obj=tmp_object, result_dir=dirname(split_object_path))
