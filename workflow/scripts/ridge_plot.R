
#### load libraries & utility function 
library("Seurat")
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
object_path <- snakemake@input[["norm_object"]]

# outputs
result_dir <- snakemake@output[["ridge_plots"]]

# parameters
step <- snakemake@wildcards[["step"]]

ab_flag <- snakemake@params[["ab_flag"]]
crispr_flag <- snakemake@params[["crispr_flag"]]
custom_flag <- snakemake@params[["custom_flag"]]

vis_gene_lists <- snakemake@config[["vis_gene_lists"]]
feature_list <- snakemake@wildcards[["feature_list"]]
category <- snakemake@wildcards[["category"]]


# ridge plot configs
width <- 4 # width of each subplot
heigth_row <- 0.5 # height of each row in the ridge plot
same.y.lims <- TRUE
assay <- "SCT"

print(category)

### load data
data_object <- readRDS(file = file.path(object_path))
Idents(object = data_object) <- category
Idents(object = data_object) <- factor(x = Idents(data_object), levels = sort(levels(data_object)))

# check if categorical metadata is all NA (can happen on data subset that does not contain the visualization category)
if(all(is.na(levels(data_object)))){
    ggsave_new(filename=paste0(cat,"_is_NA"),
               results_path=result_dir, 
               plot=ggplot() + theme_void(), 
               width=width, 
               height=height
              )
    return()   
}

# plot specs
height <- length(levels(data_object))*heigth_row + 1

# set assay, slot and fetures according to processing step and feature lists
if (step=="CORRECTED"){
    slot <- "scale.data"
    # load plotting gene list
    features <- scan(file.path(vis_gene_lists[feature_list]), character())
}else{
    slot <- "data"
    
    if(feature_list=="Antibody_Capture" | feature_list=="CRISPR_Guide_Capture" | feature_list=="Custom"){
        assay <- snakemake@config[["modality_flags"]][[feature_list]]
        features <- snakemake@config[["vis_features"]][[feature_list]]
    }else if(feature_list=="Metadata"){
        same.y.lims <- FALSE
        features <- snakemake@config[["vis_features"]][[feature_list]]
    }else{
        # load plotting gene list
        features <- scan(file.path(vis_gene_lists[feature_list]), character())
    }
}

print(feature_list)
# extract all features from the object
if(feature_list=="Metadata"){
    data_features <- colnames(data_object[[]])
}else{
    data_features <- rownames(GetAssayData(data_object, slot = slot, assay = assay))
}

# get all features of modality if indicated
if (length(features)==1){
    if(features=="all"){
        features <- data_features
    }
}

# plot only features present in the data
features <- intersect(features, data_features)

# make violin plots
for (feature in features){
    tmp_plot <- RidgePlot(object = data_object,
                                  features = feature,
                                  cols = NULL,
                                  idents = NULL,
                                  sort = FALSE,
                                  assay = assay,
                                  group.by = NULL,
                                  y.max = NULL,
                                  same.y.lims = same.y.lims,
                                  log = FALSE,
                                  ncol = 1,
                                  slot = slot,
                                  stack = FALSE,
                                  combine = TRUE,
                                  fill.by = "feature"
                                 ) + theme(legend.position = 'none')

    # save plot
    ggsave_new(filename=feature,
               results_path=result_dir, 
               plot=tmp_plot, 
               width=width, 
               height=height
              )
}
