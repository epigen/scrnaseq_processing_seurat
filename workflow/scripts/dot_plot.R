
#### load libraries & utility function 
library(Seurat)
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
object_path <- snakemake@input[["norm_object"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/merged/counts/NORMALIZED_object.rds"

# parameters
step <- "NORMALIZED"

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

# dot plot configs

height_row <- 0.5 # height of each row in the dot plot
width_col <- 0.5 # width of each col in the dot plot
result_dir <- file.path(dirname(object_path),'plots','DotPlot')

### load data
data_object <- readRDS(file = file.path(object_path))

# load plotting gene lists
gene_lists <- list()
for (gene_list_name in names(vis_gene_lists)){
    gene_lists[[gene_list_name]] <- scan(file.path(vis_gene_lists[gene_list_name]), character())
}

slot <- "data"

# make dot plots
for (cat in vis_categories){
    print(cat)
    
    Idents(object = data_object) <- cat
    Idents(object = data_object) <- factor(x = Idents(data_object), levels = sort(levels(data_object)))
    
    # check if metadata is all NA (can happen on data subset that does not contain the visualization category)
    skip <- all(is.na(levels(data_object)))
    if(!skip){
        skip <- any(is.na(Idents(object = data_object)))
    }
    
    if (any(levels(data_object)=="")){
        cells.use <- WhichCells(data_object, idents = '')
        data_object <- SetIdent(data_object, cells = cells.use, value = 'Unknown')
    }
    
    # plot specs
    height <- length(levels(data_object))*height_row + 1
    # make sure that at least the legend fits
    if (height<5){
        height <- 5
    }
    
    # plot RNA normalized/corrected expression data
    for (gene_list_name in names(gene_lists)){
        features <- unique(gene_lists[[gene_list_name]])
        
        # plot specs
        width <- width_col*length(features) + 3
        
        # plot
        if(skip){
            tmp_plot <- ggplot() + theme_void()
        }else{
            tmp_plot <- DotPlot(object = data_object,
                                  assay = "SCT",
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
                                )+ theme(axis.text.x = element_text(angle = 45, hjust=1))
        }

        
        # save plot
        ggsave_new(filename=paste0(cat,"_",gene_list_name), 
                   results_path=result_dir, 
                   plot=tmp_plot, 
                   width=width, 
                   height=height
                  )
    }
    
    # plot other modalities
    for (flag in c(ab_flag, crispr_flag, custom_flag)){
        # check if modality is used
        if(flag==''){
            next
        }
        
        # get features
        if(modality_features[[flag]][1]=='all'){
            features <- rownames(GetAssayData(data_object, slot = slot, assay = flag))
        }else{
            features <- unique(modality_features[[flag]])
        }
            
        # plot specs
        width <- width_col*length(features) + 3
        
        # plot
        if(skip){
            tmp_plot <- ggplot() + theme_void()
        }else{
            tmp_plot <- DotPlot(object = data_object,
                                  assay = flag,
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
                                )+ theme(axis.text.x = element_text(angle = 45, hjust=1))
        }
        
        # save plot
        ggsave_new(filename=paste0(cat,"_",flag), 
                   results_path=result_dir, 
                   plot=tmp_plot, 
                   width=width, 
                   height=height
                  )
    }
}