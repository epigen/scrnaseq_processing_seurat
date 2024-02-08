#### load libraries
library("Seurat")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

#### configs

# inputs
sample_object_paths <- snakemake@input #c(file.path("/nobackup/lab_bock/projects/macroIC/results/AKsmall/AK_A_transcriptome/counts/RAW_object.rds"), file.path("/nobackup/lab_bock/projects/macroIC/results/AKsmall/AK_B_transcriptome/counts/RAW_object.rds"), file.path("/nobackup/lab_bock/projects/macroIC/results/AKsmall/AK_C_transcriptome/counts/RAW_object.rds"))

# outputs
merged_object <- snakemake@output[["merged_object"]]# "/nobackup/lab_bock/projects/macroIC/results/AKsmall/merged/counts/RAW_object.rds"

# parameters
project_name <- snakemake@params[["project_name"]] #"test" 
result_dir <- dirname(merged_object)
extra_metadata_path <- snakemake@params[["extra_metadata"]]

# remove last entry if extra metadata is provided
if (extra_metadata_path != ""){
    sample_object_paths <- sample_object_paths[-length(sample_object_paths)]
}

# 'flags' for modalities
ab_flag <- snakemake@params[["ab_flag"]]#'AB'
crispr_flag <- snakemake@params[["crispr_flag"]]#'gRNA'
custom_flag <- snakemake@params[["custom_flag"]]#'HTO'


# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

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
if (extra_metadata_path != ""){
    print("extra metadata is added")
#     extra_metadata <- read.csv(extra_metadata_path, row.names = 1, header= TRUE)
    extra_metadata <- data.frame(fread(file.path(extra_metadata_path), header=TRUE), row.names=1)
    
    for (col in colnames(extra_metadata)){
        metadata_tmp <- data.frame(matrix(nrow=ncol(merged_data), ncol=1, dimnames=list(colnames(merged_data), c(col))))
        metadata_tmp[rownames(extra_metadata),col] <- extra_metadata[rownames(extra_metadata),col]
        
        # check if categorical (heuristic: less than 50 unique values) and make factor
        if (length(unique(metadata_tmp[[col]]))<50){
            metadata_tmp[, col] <- as.factor(metadata_tmp[, col])
        }
        
        merged_data[[col]] <- metadata_tmp[colnames(merged_data), col]
    }
}


### save merged data
save_seurat_object(seurat_obj=merged_data,
                   result_dir=result_dir,
                   prefix='RAW_'
                  )
