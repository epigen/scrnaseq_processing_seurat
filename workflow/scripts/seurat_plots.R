
#### load libraries & utility function 
library("Seurat")
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# inputs
object_path <- snakemake@input[["norm_object"]]

# outputs
result_dir <- snakemake@output[["plot_dir"]]

# parameters
step <- snakemake@wildcards[["step"]]
plot_type <- snakemake@wildcards[["plot_type"]]

ab_flag <- snakemake@params[["ab_flag"]]
crispr_flag <- snakemake@params[["crispr_flag"]]
custom_flag <- snakemake@params[["custom_flag"]]

vis_gene_lists <- snakemake@params[["vis_gene_lists"]]
feature_list <- snakemake@wildcards[["feature_list"]]
category <- snakemake@wildcards[["category"]]

# plot configs
if (plot_type=="VlnPlot"){
    # violin 
    pt.size <- 0 #0.1
    height <- 4 # height of each subplot
    width_col <- 0.5 # width of each col in the violin plot
    same.y.lims <- TRUE
    sorting_flag <- FALSE
} else if (plot_type=="RidgePlot"){
    # ridge
    width <- 4 # width of each subplot
    heigth_row <- 0.5 # height of each row in the ridge plot
    same.y.lims <- TRUE
    sorting_flag <- TRUE
} else if (plot_type=="DotPlot"){
    height_row <- 0.5 # height of each row in the dot plot
    width_col <- 0.5 # width of each col in the dot plot
    sorting_flag <- TRUE
} else if (plot_type=="Heatmap"){
    height_row <- 0.1 # height of each row in the heatmap plot
    width_col <- 0.001 # width of each col in the heatmap plot
    sorting_flag <- FALSE
}

# general configs (default assay)
assay <- "SCT"

print(category)

### load data
data_object <- readRDS(file = file.path(object_path))
Idents(object = data_object) <- category
Idents(object = data_object) <- factor(x = Idents(data_object), levels = sort(levels(data_object), decreasing = sorting_flag))

# check if categorical metadata is all NA (can happen on data subset that does not contain the visualization category)
if(all(is.na(levels(data_object))) | ((plot_type=="DotPlot" | plot_type=="Heatmap")&(any(is.na(Idents(object = data_object)))))){
    ggsave_new(filename=paste0(category,"_is_NA"),
               results_path=result_dir, 
               plot=ggplot() + theme_void(), 
               width=1, 
               height=1
              )
    print("Metadata is all NA.")
    quit(save = "no", status = 0)
}

# handle empy levels in DotPlots and Heatmaps
if ((plot_type=="DotPlot"  | plot_type=="Heatmap")&(any(levels(data_object)==""))){
    cells.use <- WhichCells(data_object, idents = '')
    data_object <- SetIdent(data_object, cells = cells.use, value = 'unknown')
}

# dynamic plot specs
if (plot_type=="VlnPlot"){
    width <- length(levels(data_object))*width_col + 1
} else if (plot_type=="RidgePlot"){
    height <- length(levels(data_object))*heigth_row + 1
} else if (plot_type=="DotPlot"){
    height <- max(5, length(levels(data_object))*height_row + 1)
} else if (plot_type=="Heatmap"){
    width <- length(colnames(data_object))*width_col + 2
}

# set assay, slot and fetures according to processing step and feature lists
if (step=="CORRECTED"){
    slot <- "scale.data"
    
    if (feature_list=="HVG"){
        HVG_df <- HVFInfo(object = data_object[["SCT"]], selection.method = "sct")
        features <- rownames(HVG_df[order(-HVG_df$residual_variance),])[1:100]
    }else{
        # load plotting gene list
        features <- scan(file.path(vis_gene_lists[feature_list]), character())
    }
}else{
    slot <- "data"
    
    if(feature_list=="Antibody_Capture" | feature_list=="CRISPR_Guide_Capture" | feature_list=="Custom"){
        assay <- snakemake@config[["modality_flags"]][[feature_list]]
        features <- snakemake@config[["vis_features"]][[feature_list]]
    }else if(feature_list=="Metadata"){
        same.y.lims <- FALSE
        features <- snakemake@config[["vis_features"]][[feature_list]]
    }else{
        if (feature_list=="HVG"){
            HVG_df <- HVFInfo(object = data_object[["SCT"]], selection.method = "sct")
            features <- rownames(HVG_df[order(-HVG_df$residual_variance),])[1:100]
        }else{
            # load plotting gene list
            features <- scan(file.path(vis_gene_lists[feature_list]), character())
        }
        
        # Gene expression Heatmaps always on scaled data
        if (plot_type=="Heatmap"){
            slot <- "scale.data"
        }
    }
}

print(feature_list)
# extract all features from the object
if(feature_list=="Metadata"){
    # select only numeric columns
    data_features <- colnames(data_object[[]])[sapply(data_object[[]], is.numeric)]
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
# sort features
features <- sort(features)

# check if there are any features to plot
if(length(features)==0){
    ggsave_new(filename=paste0(feature_list,"_NoFeatures"),
               results_path=result_dir, 
               plot=ggplot() + theme_void(), 
               width=1, 
               height=1
              )
    print("No features provided or present in the data.")
    quit(save = "no", status = 0)
}

# workaround for Metadata Heatmap from here: https://github.com/satijalab/seurat/issues/3366
if(plot_type=="Heatmap" & feature_list=="Metadata"){
    data_object[['metadata']] <- CreateAssayObject(data = t(x = FetchData(object = data_object, vars = features)))
    assay <- 'metadata'
    # reformat feature names due to R: “Feature names cannot have underscores ('_'), replacing with dashes ('-')”
    features <- gsub("_", "-", features)
}

# make and save Seurat plots
if (plot_type=="Heatmap"){
    # Heatmap specs
    height <- max(5, height_row*length(features) + 1)
    
    ### downsample for efficient hierarchical clustering and plotting of the same number of cells per group (ident)
    if (ncol(data_object)>30000){
        # captures count of cells for the ident with the fewest (but at least 100)
        maxcells  <- min(100,table(Idents(data_object)))
        # downsample
        set.seed(42)
        data_object <- subset(data_object, downsample = maxcells)
        # update heatmap width accordingly
        width <- length(colnames(data_object))*width_col + 4
    }
    
    # cluster features (rows) using hclust
    tmp_data <- as.data.frame(GetAssayData(object = data_object, slot = slot, assay = assay))
    # subset by features
    tmp_data <- tmp_data[features,]
    # cluster using hclust
    hc <- hclust(dist(tmp_data))
    
    # make Heatmap
    tmp_plot <- DoHeatmap(object = data_object,
                                  features = features[hc$order],# reorder features by hclust result
                                  cells = NULL,
                                  group.by = "ident",
                                  group.bar = TRUE,
                                  group.colors = NULL,
                                  disp.min = -2.5,
                                  disp.max = NULL,
                                  slot = slot,
                                  assay = assay,
                                  label = TRUE,
                                  size = 2,
                                  hjust = 0.5,
                                  angle = 0,
                                  raster = TRUE,
                                  draw.lines = TRUE,
                                  lines.width = 20,
                                  group.bar.height = 0.02,
                                  combine = TRUE
                                ) + guides(colour="none") + 
    scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick2", midpoint = 0, space ="Lab", na.value = "white")
#     scale_fill_gradientn(colors = c("royalblue4", "white", "firebrick2"), na.value = "white")

    # save plot
    ggsave_new(filename=feature_list,
               results_path=result_dir, 
               plot=tmp_plot, 
               width=width, 
               height=height
              )
} else if (plot_type=="DotPlot"){
    # DotPlot specs
    width <- width_col*length(features) + 3
    # make DotPlot
    tmp_plot <- DotPlot(object = data_object,
                                  assay = assay,
                                  features = features,
                                  cols = c("lightgrey", "blue"),
                                  col.min = -2.5,
                                  col.max = 2.5,
                                  dot.min = 0,
                                  dot.scale = 6,
                                  idents = NULL,
                                  group.by = NULL,
                                  split.by = NULL,
                                  cluster.idents = FALSE,
                                  scale = TRUE,
                                  scale.by = "radius",
                                  scale.min = NA,
                                  scale.max = NA
                                ) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    scale_color_gradient2(midpoint=0, low="royalblue4", mid="grey80", high="firebrick2", space ="Lab")
    # save plot
    ggsave_new(filename=feature_list,
               results_path=result_dir, 
               plot=tmp_plot, 
               width=width, 
               height=height
              )
} else if (plot_type=="VlnPlot"){
    for (feature in features){
        # make violin plots
        tmp_plot <- VlnPlot(object = data_object,
                        features = feature,
                        cols = NULL,
                        pt.size = pt.size,
                        idents = NULL,
                        sort = FALSE,
                        assay = assay,
                        group.by = NULL,
                        split.by = NULL,
                        adjust = 1,
                        y.max = NULL,
                        same.y.lims = same.y.lims,
                        log = FALSE,
                        ncol = 1,#ncols,
                        slot = slot,
                        split.plot = FALSE,
                        stack = FALSE,
                        combine = TRUE,
                        fill.by = "feature",
                        flip = FALSE
                        ) + theme(legend.position = 'none')
        
        # save plot
        ggsave_new(filename=feature,
                   results_path=result_dir, 
                   plot=tmp_plot, 
                   width=width, 
                   height=height
                  )
    }
} else if (plot_type=="RidgePlot"){
    for (feature in features){
        # make ridge plots
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
}
