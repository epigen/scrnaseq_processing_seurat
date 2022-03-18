
#### load libraries & utility function 
library(Seurat)
# source utility functions
source("workflow/scripts/utils.R")

#### configs

# inputs
sample_dir <- snakemake@input[[1]]#file.path("/research/lab_bock/projects/macroIC/data/BSA_0560_MF_AK_MPH/OUT/COUNT/AK_A_transcriptome/")
metadata_path <- snakemake@input[[2]]#file.path("/research/lab_bock/projects/macroIC/data/BSA_0560_MF_AK_MPH/OUT/COUNT/AK_A_transcriptome/", "QC_categories.csv")

# outputs
result_object <- snakemake@output[["sample_object"]]# file.path("/nobackup/lab_bock/projects/macroIC/results/AK_A_transcriptome/counts/RAW_object.rds")

# parameters
saveCounts <- snakemake@params[["saveCounts"]]['raw']

sample_name <- snakemake@params[["sample"]]#"AK_A_transcriptome"
result_dir <- dirname(result_object)
# 'flags' for modalities
ab_flag <- snakemake@params[["ab_flag"]]#'AB'
crispr_flag <- snakemake@params[["crispr_flag"]]#'gRNA'
custom_flag <- snakemake@params[["custom_flag"]]#'HTO'
# gRNA assignment threshold
crispr_umi_threshold <- snakemake@params[["crispr_umi_threshold"]]#1
# gRNA to gene REGEX
grna_regex = snakemake@params[["grna_regex"]]#"\\-."
# REGEX-based metadata extension for Seurat::PercentageFeatureSet()
percentage_regex = snakemake@params[["percentage_regex"]]#list(mito.percent="^MT-")
# eval based metadata column transformations
metadata_eval = snakemake@params[["metadata_eval"]]


# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

#### load data (future: check here for 10X/count matrix flag)
data <- Read10X(data.dir = file.path(sample_dir, "filtered_feature_bc_matrix"))

# Initialize the Seurat object with the raw data
seurat_obj <- CreateSeuratObject(counts = data$'Gene Expression', project=sample_name)

if(ab_flag!=''){
    # create a new assay to store Antibody information
    seurat_obj[[ab_flag]] <- CreateAssayObject(counts = data$'Antibody Capture')
}

if(crispr_flag!=''){
    # create a new assay to store CRISPR guide information
    seurat_obj[[crispr_flag]] <- CreateAssayObject(counts = data$'CRISPR Guide Capture')
}

if(custom_flag!=''){
    # create a new assay to store Custom information
    seurat_obj[[custom_flag]] <- CreateAssayObject(counts = data$'Custom')
}


#### load & add metadata
metadata <- read.csv(metadata_path, row.names = 1, header= TRUE)
metadata$batch <- sample_name

for (col in colnames(metadata)){
    seurat_obj[[col]] <- metadata[colnames(seurat_obj),col]
}

# extend metadata with Seurat::PercentageFeatureSet
for (name in names(percentage_regex)){
    seurat_obj[[name]] <- PercentageFeatureSet(seurat_obj, pattern = percentage_regex[[name]])
}

if(crispr_flag!=''){
    # gRNA and KO assignment
    grna_data <- GetAssayData(object = seurat_obj, slot = "counts", assay = crispr_flag)

    # load threshold data from cellranger
    thresholds <- read.csv(file.path(sample_dir, "crispr_analysis", "protospacer_umi_thresholds.csv"), row.names = 1, header= TRUE)

    # substitute '_' with '-' (R constraint: “Feature names cannot have underscores ('_'), replacing with dashes ('-')”)
    rownames(thresholds) <- gsub("_", "-", rownames(thresholds))

    for (grna in rownames(grna_data)){
        if (grna %in% rownames(thresholds)){
            grna_data[grna, grna_data[grna,] < thresholds[grna,'UMI.threshold']] <- 0
            grna_data[grna, grna_data[grna,] <= crispr_umi_threshold] <- 0
        } else{
            grna_data[grna, ] <- 0
            sprintf("in sample %s 0 UMIs of guide RNA %s found", sample_name, grna)
        }
    }

    # perfrom gRNA and KO target assignment
    protospacer_calls_df <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(c(),c("guide_call", "KO_call"))))

    for (barcode in colnames(seurat_obj)){

        # guide RNA and target gene can be assigned
        if (sum(grna_data[,barcode]>0)==1){
            protospacer_calls_df[barcode,"guide_call"] <- rownames(grna_data)[grna_data[,barcode]>0]
            protospacer_calls_df[barcode,"KO_call"] <- gsub(grna_regex,"",protospacer_calls_df[barcode,"guide_call"])
        }

        # multiple guide RNAs have been detected
        if (sum(grna_data[,barcode]>0)>1){
            protospacer_calls_df[barcode,"guide_call"] <- "Multiplet"
            protospacer_calls_df[barcode,"KO_call"] <- "Multiplet"
        }

        # no guide RNAs have been detected
        if (sum(grna_data[,barcode]>0)==0){
            protospacer_calls_df[barcode,"guide_call"] <- "Negative"
            protospacer_calls_df[barcode,"KO_call"] <- "Negative"
        }
    }

    # add guide and KO assignment to metadata
    seurat_obj[["guide_call"]] <- protospacer_calls_df[colnames(seurat_obj), 'guide_call']
    seurat_obj[["KO_call"]] <- protospacer_calls_df[colnames(seurat_obj), 'KO_call']
}

# add eval based metadata columns from config
for (name in names(metadata_eval)){
    metadata <- seurat_obj[[]] # ensure that cells ordered correctly and current/latest metadtata is used
    seurat_obj[[name]] <- eval(parse(text=metadata_eval[name]))
}

#### SAVE RESULTS
save_seurat_object(seurat_obj=seurat_obj, result_dir=result_dir, prefix='RAW_', ab_flag=ab_flag, crispr_flag=crispr_flag, custom_flag=custom_flag, saveCounts=saveCounts)

