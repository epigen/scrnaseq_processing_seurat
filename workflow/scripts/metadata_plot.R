
#### load libraries & utility function 
library("inspectdf")
library(tibble)
library("ggplot2")
# source utility functions
source("workflow/scripts/utils.R")

# inputs
metadata_path <- snakemake@input[["metadata"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/merged/counts/NORMALIZED_metadata.csv"

# parameters
step <- snakemake@params[["step"]] #"NORMALIZED"


# load data
metadata <- read.csv(metadata_path, row.names = 1, header= TRUE)
metadata <- as_tibble(metadata)

### plot metadata types
metadata_type_stats <- inspect_types(metadata)
tmp_plot <- show_plot(metadata_type_stats)

# plots specs
width <- 10
height <- 10

# save plot
ggsave_new(filename=paste0(step,"_metadata_types"), 
           results_path=file.path(dirname(metadata_path), 'plots'), 
           plot=tmp_plot, 
           width=width, 
           height=height
          )

### categorical data
# inspect categorical data
metadata_cat_stats <- inspect_cat(metadata)

# plot stats
tmp_plot <- metadata_cat_stats %>% show_plot(label_thresh=0.01)

# plots specs
height <- nrow(metadata_cat_stats)*0.5

# save plot
ggsave_new(filename=paste0(step,"_metadata_categorical"), 
           results_path=file.path(dirname(metadata_path), 'plots'), 
           plot=tmp_plot, 
           width=width, 
           height=height
          )


### numerical data
# inspect numerical data
metadata_num_stats <- inspect_num(metadata, breaks=100)

# plot stats
tmp_plot <- metadata_num_stats %>% show_plot()

# plots specs
height <- nrow(metadata_num_stats)/3*1.5

# save plot
ggsave_new(filename=paste0(step,"_metadata_numerical"), 
           results_path=file.path(dirname(metadata_path), 'plots'), 
           plot=tmp_plot, 
           width=width, 
           height=height
          )


### save all statistics as CSV files
# make result directory if not exist
result_dir <- file.path(dirname(metadata_path), 'stats')
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

for (cat in names(metadata_cat_stats$levels)){
    write.csv(metadata_cat_stats$levels[[cat]], file=file.path(result_dir, paste0(step,"_metadata_",cat,".csv")), row.names=FALSE)
}

write.csv(metadata_num_stats[,-ncol(metadata_num_stats)], file=file.path(result_dir, paste0(step,"_metadata_","numerical",".csv")), row.names=FALSE)

# library(openxlsx)
# write.xlsx(metadata_cat_stats$levels, file = file.path(result_dir, paste0(step,"_metadata_stats",".xlsx")))