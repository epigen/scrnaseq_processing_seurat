
#### load libraries & utility function 
library("Seurat")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# helper function to assign each cell's gRNA and KO calls
assign_grna_KO <- function(col) {
    non_zero_indices <- which(col != 0)
    non_zero_guides <- names(non_zero_indices)
    non_zero_kos <- unique(gsub(grna_regex,"",non_zero_guides))
    
    if (length(non_zero_guides) > 0) {
      gRNAcall <- paste(non_zero_guides, collapse = "_")
      KOcall <- paste(non_zero_kos, collapse = "_")
  } else {
      # in case there are no non-zero values
      gRNAcall <- NA 
      KOcall <- NA
  }
    gRNA_n <- length(non_zero_guides)
    KO_n <- length(non_zero_kos)
  
  return(c(gRNAcall = gRNAcall, gRNA_n = gRNA_n, KOcall = KOcall, KO_n = KO_n))
}

#### configs

# inputs
sample_dir <- snakemake@input[[1]]
metadata_path <- snakemake@params[["metadata"]]

# outputs
result_object <- snakemake@output[["sample_object"]]

# parameters
sample_name <- snakemake@wildcards[["sample"]]
result_dir <- dirname(result_object)
# 'flags' for modalities
ab_flag <- snakemake@config[["modality_flags"]][['Antibody_Capture']]
crispr_flag <- snakemake@config[["modality_flags"]][['CRISPR_Guide_Capture']]
custom_flag <- snakemake@config[["modality_flags"]][['Custom']]
# gRNA assignment threshold
crispr_umi_threshold <- snakemake@config[["crispr_umi_threshold"]]
# gRNA to gene REGEX
grna_regex <- snakemake@config[["grna_regex"]]
# REGEX-based metadata extension for Seurat::PercentageFeatureSet()
percentage_regex <- snakemake@config[["percentage_regex"]]
# eval based metadata column transformations
metadata_eval <- snakemake@config[["metadata_eval"]]


#### load data & Initialize the Seurat object with the raw data
if (dir.exists(file.path(sample_dir, "filtered_feature_bc_matrix"))){
    print("Load 10X data")
    is10x <- TRUE
    data <- Read10X(data.dir = file.path(sample_dir, "filtered_feature_bc_matrix"))
    seurat_obj <- CreateSeuratObject(counts = data$'Gene Expression', project=sample_name)
}else{
    print("Load mtx data")
    is10x <- FALSE
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
    # read the features file again to get our modality 'type' column
    features <- read.table(
        file.path(sample_dir, "features.tsv"),
        sep = "\t",
        header = FALSE,
        col.names = c("id", "name", "type")
    )
    # The features read by ReadMtx become the rownames, so we use them to align
    rownames(features) <- rownames(data)
    rna_features <- features[features$type == "Gene Expression", ]
    rna_counts <- data[rownames(rna_features), ]
    
    seurat_obj <- CreateSeuratObject(counts = rna_counts, project = sample_name)
}

print(seurat_obj)

if(is10x){
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
}else{
    if(ab_flag!=''){
        # create a new assay to store Antibody information
        ab_features <- features[features$type == "Antibody Capture", ]
        ab_counts <- data[rownames(ab_features), , drop = FALSE]
        seurat_obj[[ab_flag]] <- CreateAssayObject(counts = ab_counts)
    }
    
    if(crispr_flag!=''){
        # create a new assay to store CRISPR guide information
        gdo_features <- features[features$type == "CRISPR Guide Capture", ]
        gdo_counts <- data[rownames(gdo_features), , drop = FALSE]
        seurat_obj[[crispr_flag]] <- CreateAssayObject(counts = gdo_counts)
    }
    
    if(custom_flag!=''){
        # create a new assay to store Custom information
        custom_features <- features[features$type == "Custom", ]
        custom_counts <- data[rownames(custom_features), , drop = FALSE]
        seurat_obj[[custom_flag]] <- CreateAssayObject(counts = custom_counts)
    }
}

#### load & add metadata (if exists)
if (length(metadata_path) > 0 ){
    if (file.exists(file.path(metadata_path))){
        print("load & add metadata")
        
#         metadata <- read.csv(file.path(metadata_path), row.names = 1, header= TRUE)
        metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)

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

# assign gRNA and KO calls
if(crispr_flag!=''){
    # gRNA and KO assignment
    grna_data <- GetAssayData(object = seurat_obj, slot = "counts", assay = crispr_flag)

    if(is10x){
        # load threshold data from cellranger output
        thresholds <- data.frame(fread(file.path(sample_dir, "crispr_analysis", "protospacer_umi_thresholds.csv"), header=TRUE), row.names=1)
    
        # substitute '_' with '-' (R constraint: “Feature names cannot have underscores ('_'), replacing with dashes ('-')”)
        rownames(thresholds) <- gsub("_", "-", rownames(thresholds))
    
        # apply cellranger and configured thresholds
        for (grna in rownames(grna_data)){
            
            if (grna %in% rownames(thresholds)){
                grna_data[grna, grna_data[grna,] < thresholds[grna,'UMI.threshold']] <- 0
            } else{
                grna_data[grna, ] <- 0
                sprintf("in sample %s no UMI threshold for guide RNA %s was found -> setting it to zero", sample_name, grna)
            }
            # apply configured UMI threshold
            grna_data[grna, grna_data[grna,] <= crispr_umi_threshold] <- 0
        }
    }else{
        # apply only the configured thresholds
        for (grna in rownames(grna_data)){
             grna_data[grna, grna_data[grna,] <= crispr_umi_threshold] <- 0
        }
    }

    # perfrom gRNA and KO target assignment
    protospacer_calls_df <- as.data.frame(t(apply(grna_data, 2, assign_grna_KO)))
    protospacer_calls_df$gRNA_n <- as.numeric(protospacer_calls_df$gRNA_n)
    protospacer_calls_df$KO_n <- as.numeric(protospacer_calls_df$KO_n)
    
#     protospacer_calls_df <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(c(),c("gRNAcall", "KOcall"))))

#     for (barcode in colnames(seurat_obj)){

#         # guide RNA and target gene can be assigned
#         if (sum(grna_data[,barcode]>0)==1){
#             protospacer_calls_df[barcode,"gRNAcall"] <- rownames(grna_data)[grna_data[,barcode]>0]
#             protospacer_calls_df[barcode,"KOcall"] <- gsub(grna_regex,"",protospacer_calls_df[barcode,"gRNAcall"])
#         }

#         # multiple guide RNAs have been detected
#         if (sum(grna_data[,barcode]>0)>1){
#             protospacer_calls_df[barcode,"gRNAcall"] <- "Multiplet"
#             protospacer_calls_df[barcode,"KOcall"] <- "Multiplet"
#         }

#         # no guide RNAs have been detected
#         if (sum(grna_data[,barcode]>0)==0){
#             protospacer_calls_df[barcode,"gRNAcall"] <- "Negative"
#             protospacer_calls_df[barcode,"KOcall"] <- "Negative"
#         }
#     }

    # add guide and KO assignment to metadata
    for (col in colnames(protospacer_calls_df)){
        seurat_obj[[col]] <- protospacer_calls_df[colnames(seurat_obj), col]
    }
}

# add eval based metadata columns from config
for (name in names(metadata_eval)){
    metadata <- seurat_obj[[]] # ensure that cells ordered correctly and current/latest metadtata is used
    seurat_obj[[name]] <- eval(parse(text=metadata_eval[name]))
}

#### SAVE RESULTS
save_seurat_object(seurat_obj=seurat_obj,
                   result_dir=result_dir
                  )

