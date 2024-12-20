
#### load libraries & utility function 
library("inspectdf")
library("tibble")
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# inputs
metadata_path <- snakemake@input[["metadata"]]

# outputs
result_dir <- snakemake@output[["metadata_plots"]]
stats_dir <- snakemake@output[["metadata_stats"]]

# load data
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)
metadata <- as_tibble(metadata)

### plot metadata types
metadata_type_stats <- inspect_types(metadata)
tmp_plot <- show_plot(metadata_type_stats)

# plots specs
width <- 5
height <- 5

# save plot
ggsave_new(filename="types", 
           results_path=result_dir, 
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
width <- 10

# save plot
ggsave_new(filename="categorical", 
           results_path=result_dir, 
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
ggsave_new(filename="numerical",
           results_path=result_dir,
           plot=tmp_plot, 
           width=width, 
           height=height
          )


### save all statistics as CSV files
dir.create(stats_dir, recursive = TRUE)

for (cat in names(metadata_cat_stats$levels)){
#     write.csv(metadata_cat_stats$levels[[cat]], file=file.path(stats_dir, paste0(step,"_metadata_",cat,".csv")), row.names=FALSE)
    fwrite(as.data.frame(metadata_cat_stats$levels[[cat]]), file=file.path(stats_dir, paste0("metadata_",cat,".csv")), row.names=FALSE)
}

# write.csv(metadata_num_stats[,-ncol(metadata_num_stats)], file=file.path(stats_dir, paste0(step,"_metadata_","numerical",".csv")), row.names=FALSE)
fwrite(as.data.frame(metadata_num_stats[,-ncol(metadata_num_stats)]), file=file.path(stats_dir, "metadata_numerical.csv"), row.names=FALSE)
