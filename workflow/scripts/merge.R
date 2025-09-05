#### load libraries
library("Seurat")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

#### configs

# inputs
sample_object_paths <- snakemake@input[["samples"]]
extra_metadata_path <- snakemake@input[["extra_metadata"]]

# outputs
merged_object <- snakemake@output[["merged_object"]]

# parameters
project_name <- snakemake@config[["project_name"]] #"test" 
result_dir <- dirname(merged_object)

# 'flags' for modalities
ab_flag <- snakemake@config[["modality_flags"]][['Antibody_Capture']]
crispr_flag <- snakemake@config[["modality_flags"]][['CRISPR_Guide_Capture']]
custom_flag <- snakemake@config[["modality_flags"]][['Custom']]


### load data
sample_objects = c()
for (sample_path in sample_object_paths){
        sample_objects <- append(sample_objects, readRDS(file = file.path(sample_path)))
}

# if more than one sample merge into one object, otherwise just rename
if (length(sample_objects)>1){
    merged_data <- merge(sample_objects[[1]],
                         y = sample_objects[2:length(sample_objects)],
                         project = project_name,
                         add.cell.ids = unlist(lapply(sample_objects, Project))
              )
}else{
    merged_data <- sample_objects[[1]]
}

# if extra metadata is provided add it to the merged object
if (length(extra_metadata_path)!=0){
    print("extra metadata is added")
    extra_metadata <- data.frame(fread(file.path(extra_metadata_path), header=TRUE), row.names=1)
    
    for (col in colnames(extra_metadata)){
        metadata_tmp <- data.frame(matrix(nrow=ncol(merged_data), ncol=1, dimnames=list(colnames(merged_data), c(col))))
        metadata_tmp[rownames(extra_metadata),col] <- extra_metadata[rownames(extra_metadata),col]
        
        # check if categorical (heuristic: less than 50 unique values) and make factor
        if (length(unique(metadata_tmp[[col]]))<50){
            metadata_tmp[, col] <- as.factor(metadata_tmp[, col])
            
            # check if any entry is an empty string i.e., "" and replace with "unknown"
            metadata_tmp[[col]][metadata_tmp[[col]] == ""] <- "unknown"
        }
        
        merged_data[[col]] <- metadata_tmp[colnames(merged_data), col]
    }
}

### save merged data
save_seurat_object(seurat_obj=merged_data,result_dir=dirname(merged_object))
