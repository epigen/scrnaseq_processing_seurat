#### load libraries
library(Seurat)
# source utility functions
source("workflow/scripts/utils.R")

#### configs

# inputs
sample_object_paths <- snakemake@input #c(file.path("/nobackup/lab_bock/projects/macroIC/results/AKsmall/AK_A_transcriptome/counts/RAW_object.rds"), file.path("/nobackup/lab_bock/projects/macroIC/results/AKsmall/AK_B_transcriptome/counts/RAW_object.rds"), file.path("/nobackup/lab_bock/projects/macroIC/results/AKsmall/AK_C_transcriptome/counts/RAW_object.rds"))

# outputs
merged_object <- snakemake@output[["merged_object"]]# "/nobackup/lab_bock/projects/macroIC/results/AKsmall/merged/counts/RAW_object.rds"

# parameters
saveCounts <- snakemake@params[["saveCounts"]]['raw']

project_name <- snakemake@params[["project_name"]] #"test" 
result_dir <- dirname(merged_object)

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

### save merged data
save_seurat_object(seurat_obj=merged_data, result_dir=result_dir, prefix='RAW_', ab_flag=ab_flag, crispr_flag=crispr_flag, custom_flag=custom_flag, saveCounts=saveCounts)
