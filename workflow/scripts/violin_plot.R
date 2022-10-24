
#### load libraries & utility function 
library(Seurat)
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
object_path <- snakemake@input[["norm_object"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/merged/counts/NORMALIZED_object.rds"

# parameters
step <- snakemake@params[["step"]] #"NORMALIZED"

ab_flag <- snakemake@params[["ab_flag"]]#'AB'
crispr_flag <- snakemake@params[["crispr_flag"]]#'gRNA'
custom_flag <- snakemake@params[["custom_flag"]]#'HTO'

vis_categories <- snakemake@params[["vis_categories"]] #c('condition','KO_call')#,'condKO')
vis_gene_lists <- snakemake@params[["vis_gene_lists"]] #list(IC_genes = "/research/home/sreichl/projects/macroIC/metadata/marker_genes/IC_genes.txt", immune_genes = "/research/home/sreichl/projects/macroIC/metadata/marker_genes/immune_genes.txt")

modality_features <- list()
modality_features[[ab_flag]] <- snakemake@params[["ab_features"]] #c('all')
modality_features[[crispr_flag]] <- snakemake@params[["crispr_features"]] #c('all')
modality_features[[custom_flag]] <- snakemake@params[["custom_features"]] #c('HTO-THP1-A-untreated','HTO-THP1-A-8h-cytokines','HTO-THP1-A-24h-cytokines') #c('all')
modality_features[["Metadata"]] <- snakemake@params[["metadata_features"]]

# options(repr.plot.width=width, repr.plot.height=height)

# violin plot configs
pt.size <- 0 #0.1
ncols <- 10 #number of columns
height_plot <- 4 # height of each subplot
width_col <- 0.5 # width of each col in the violin plot
result_dir <- file.path(dirname(object_path),'plots')

### load data
data_object <- readRDS(file = file.path(object_path))

# load plotting gene lists
gene_lists <- list()
for (gene_list_name in names(vis_gene_lists)){
    gene_lists[[gene_list_name]] <- scan(file.path(vis_gene_lists[gene_list_name]), character())
}

if (step=="NORMALIZED"){
    slot <- "data"
}else{
    slot <- "scale.data"
}

# make violin plots
for (cat in vis_categories){
    print(cat)
    Idents(object = data_object) <- cat
    Idents(object = data_object) <- factor(x = Idents(data_object), levels = sort(levels(data_object)))
    
    # check if metadata is all NA (can happen on data subset that does not contain the visualization category)
    skip <- all(is.na(levels(data_object)))
    
    # plot specs
    width_plot <- length(levels(data_object))*width_col + 1
    width <- width_plot*ncols
    
    # plot RNA normalized/corrected expression data
    for (gene_list_name in names(gene_lists)){
        print(gene_list_name)
        features <- gene_lists[[gene_list_name]]
        
        # plot specs
        height <- height_plot*ceiling(length(features)/ncols)
        
        # plot
        if(skip){
            tmp_plot <- ggplot() + theme_void()
        }else{
            tmp_plot <- VlnPlot(object = data_object,
                            features = features,
                            cols = NULL,
                            pt.size = pt.size,
                            idents = NULL,
                            sort = FALSE,
                            assay = "SCT",
                            group.by = NULL,
                            split.by = NULL,
                            adjust = 1,
                            y.max = NULL,
                            same.y.lims = TRUE,
                            log = FALSE,
                            ncol = ncols,
                            slot = slot,
                            split.plot = FALSE,
                            stack = FALSE,
                            combine = TRUE,
                            fill.by = "feature",
                            flip = FALSE
                            )
        }
        
        # save plot
        ggsave_new(filename=paste0(step,"_violin_plot_",cat,"_",gene_list_name), 
                   results_path=result_dir, 
                   plot=tmp_plot, 
                   width=width, 
                   height=height
                  )
    }
    
    # other modalities are not CORRECTED ie no added value in plotting NORMALIZED data again
    if (step=="CORRECTED"){
        next
    }
    
    # plot other modalities
    for (flag in c(ab_flag, crispr_flag, custom_flag, "Metadata")){
        print(flag)
        # check if modality is used
        if(flag==''){
            next
        }
        # check if modality features are provided
        if(length(modality_features[[flag]])<1){
            next
        }
        
        # get features
        if(modality_features[[flag]][1]=='all'){
            features <- rownames(GetAssayData(data_object, slot = slot, assay = flag))
        }else{
            features <- modality_features[[flag]]
        }
        
        if(flag=='Metadata'){
            assay <- 'SCT'
            same.y.lims <- FALSE
        }else{
            assay <- flag
            same.y.lims <- TRUE
        }
            
        # plot specs
        height <- height_plot*ceiling(length(features)/ncols)
        
        # plot
        if(skip){
            tmp_plot <- ggplot() + theme_void()
        }else{
            tmp_plot <- VlnPlot(object = data_object,
                                features = features,
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
                                ncol = ncols,
                                slot = slot,
                                split.plot = FALSE,
                                stack = FALSE,
                                combine = TRUE,
                                fill.by = "feature",
                                flip = FALSE
                                )
        }
        
        # save plot
        ggsave_new(filename=paste0(step,"_violin_plot_",cat,"_",flag), 
                   results_path=result_dir, 
                   plot=tmp_plot, 
                   width=width, 
                   height=height
                  )
    }
}