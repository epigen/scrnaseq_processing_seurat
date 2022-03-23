
#### load libraries & utility function 
library(Seurat)
library("ggplot2")
# source utility functions
source("workflow/scripts/utils.R")

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

# options(repr.plot.width=width, repr.plot.height=height)

# heatmap plot configs

height_row <- 0.1 # height of each row in the heatmap plot
width_col <- 0.001 # width of each col in the heatmap plot
result_dir <- file.path(dirname(object_path),'plots')

### load data
data_object <- readRDS(file = file.path(object_path))

# load plotting gene lists
gene_lists <- list()
for (gene_list_name in names(vis_gene_lists)){
    gene_lists[[gene_list_name]] <- scan(file.path(vis_gene_lists[gene_list_name]), character())
}

HVG_df <- HVFInfo(object = data_object, assay = "SCT")
gene_lists[["HVG100"]] <- rownames(HVG_df[order(-HVG_df$residual_variance),])[1:100]

# make heatmap plots
for (cat in vis_categories){
    Idents(object = data_object) <- cat
    Idents(object = data_object) <- factor(x = Idents(data_object), levels = sort(levels(data_object)))
    
    # plot specs
    width <- length(colnames(data_object))*width_col + 2
    
    # plot RNA normalized/corrected expression data
    for (gene_list_name in names(gene_lists)){
        slot <- "scale.data"
        features <- gene_lists[[gene_list_name]]
        
        # plot specs
        height <- height_row*length(features) + 1
        # make sure that at least the legend fits
        if (height<5){
            height <- 5
        }
        
        # plot
        tmp_plot <- DoHeatmap(object = data_object,
                              features = features,
                              cells = NULL,
                              group.by = "ident",
                              group.bar = TRUE,
                              group.colors = NULL,
                              disp.min = -2.5,
                              disp.max = NULL,
                              slot = slot,
                              assay = "SCT",
                              label = TRUE,
                              size = 2,
                              hjust = 0.5,
                              angle = 0,
                              raster = TRUE,
                              draw.lines = TRUE,
                              lines.width = NULL,
                              group.bar.height = 0.02,
                              combine = TRUE
                            )
        
        # save plot
        ggsave_new(filename=paste0(step,"_heatmap_",cat,"_",gene_list_name), 
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
    for (flag in c(ab_flag, crispr_flag, custom_flag)){
        slot <- "data"
        # check if modality is used
        if(flag==''){
            next
        }
        
        # get features
        if(modality_features[[flag]][1]=='all'){
            features <- rownames(GetAssayData(data_object, slot = slot, assay = flag))
        }else{
            features <- modality_features[[flag]]
        }
            
        # plot specs
        height <- height_row*length(features) + 1
        # make sure that at least the legend fits
        if (height<5){
            height <- 5
        }
        
        
        # plot
        tmp_plot <- DoHeatmap(object = data_object,
                              features = features,
                              cells = NULL,
                              group.by = "ident",
                              group.bar = TRUE,
                              group.colors = NULL,
                              disp.min = -2.5,
                              disp.max = NULL,
                              slot = slot,
                              assay = flag,
                              label = TRUE,
                              size = 3,
                              hjust = 0.5,
                              angle = 0,
                              raster = TRUE,
                              draw.lines = TRUE,
                              lines.width = NULL,
                              group.bar.height = 0.02,
                              combine = TRUE
                            )
        
        # save plot
        ggsave_new(filename=paste0(step,"_heatmap_",cat,"_",flag), 
                   results_path=result_dir, 
                   plot=tmp_plot, 
                   width=width, 
                   height=height
                  )
    }
}