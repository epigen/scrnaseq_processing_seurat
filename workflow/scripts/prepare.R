
#### load libraries & utility function 
library(Seurat)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

#### configs

# inputs
sample_dir <- snakemake@input[[1]]#file.path("/research/lab_bock/projects/macroIC/data/BSA_0560_MF_AK_MPH/OUT/COUNT/AK_A_transcriptome/")
metadata_path <- snakemake@input[[2]]#file.path("/research/lab_bock/projects/macroIC/data/BSA_0560_MF_AK_MPH/OUT/COUNT/AK_A_transcriptome/", "QC_categories.csv")

# outputs
result_object <- snakemake@output[["sample_object"]]# file.path("/nobackup/lab_bock/projects/macroIC/results/AK_A_transcriptome/counts/RAW_object.rds")

# parameters
sample_name <- snakemake@params[["sample"]]#"AK_A_transcriptome"
result_dir <- dirname(result_object)
# 'flags' for modalities
ab_flag <- snakemake@params[["ab_flag"]]#'AB'
crispr_flag <- snakemake@params[["crispr_flag"]]#'gRNA'
custom_flag <- snakemake@params[["custom_flag"]]#'HTO'
# gRNA assignment threshold
crispr_umi_threshold <- snakemake@params[["crispr_umi_threshold"]]#1
# gRNA to gene REGEX
grna_regex <- snakemake@params[["grna_regex"]]#"\\-."
# REGEX-based metadata extension for Seurat::PercentageFeatureSet()
percentage_regex <- snakemake@params[["percentage_regex"]]#list(mito.percent="^MT-")
# eval based metadata column transformations
metadata_eval <- snakemake@params[["metadata_eval"]]


# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

#### load data & Initialize the Seurat object with the raw data
if (dir.exists(file.path(sample_dir, "filtered_feature_bc_matrix"))){
    print("Load 10X data")
    data <- Read10X(data.dir = file.path(sample_dir, "filtered_feature_bc_matrix"))
    seurat_obj <- CreateSeuratObject(counts = data$'Gene Expression', project=sample_name)
}else{
    print("Load mtx data")
    data <- ReadMtx(
        mtx = file.path(sample_dir, "matrix.mtx"),
        cells = file.path(sample_dir, "barcodes.tsv"),
        features = file.path(sample_dir, "features.tsv"),
        cell.column = 1,
        feature.column = 1,
        cell.sep = "\t",
        feature.sep = "\t",
        skip.cell = 0,
        skip.feature = 0,
        mtx.transpose = FALSE,
        unique.features = TRUE,
        strip.suffix = FALSE
    )
    seurat_obj <- CreateSeuratObject(counts = data, project=sample_name)
}

print(seurat_obj)


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


#### load & add metadata (if exists)
if (length(metadata_path) > 0 ){
    if (file.exists(file.path(metadata_path))){
        print("load & add metadata")
        
        metadata <- read.csv(file.path(metadata_path), row.names = 1, header= TRUE)

        for (col in colnames(metadata)){
            seurat_obj[[col]] <- metadata[colnames(seurat_obj),col]
        }
    }
}

# add batch info to metadata
seurat_obj$batch <- sample_name

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
    protospacer_calls_df <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(c(),c("gRNAcall", "KOcall"))))

    for (barcode in colnames(seurat_obj)){

        # guide RNA and target gene can be assigned
        if (sum(grna_data[,barcode]>0)==1){
            protospacer_calls_df[barcode,"gRNAcall"] <- rownames(grna_data)[grna_data[,barcode]>0]
            protospacer_calls_df[barcode,"KOcall"] <- gsub(grna_regex,"",protospacer_calls_df[barcode,"gRNAcall"])
        }

        # multiple guide RNAs have been detected
        if (sum(grna_data[,barcode]>0)>1){
            protospacer_calls_df[barcode,"gRNAcall"] <- "Multiplet"
            protospacer_calls_df[barcode,"KOcall"] <- "Multiplet"
        }

        # no guide RNAs have been detected
        if (sum(grna_data[,barcode]>0)==0){
            protospacer_calls_df[barcode,"gRNAcall"] <- "Negative"
            protospacer_calls_df[barcode,"KOcall"] <- "Negative"
        }
    }

    # add guide and KO assignment to metadata
    seurat_obj[["gRNAcall"]] <- protospacer_calls_df[colnames(seurat_obj), 'gRNAcall']
    seurat_obj[["KOcall"]] <- protospacer_calls_df[colnames(seurat_obj), 'KOcall']
}

# add eval based metadata columns from config
for (name in names(metadata_eval)){
    metadata <- seurat_obj[[]] # ensure that cells ordered correctly and current/latest metadtata is used
    seurat_obj[[name]] <- eval(parse(text=metadata_eval[name]))
}

#### SAVE RESULTS
save_seurat_object(seurat_obj=seurat_obj,
                   result_dir=result_dir,
                   prefix='prep_'
                  )

