### Utility Functions

# save Seurat object results
save_seurat_object <- function (seurat_obj, result_dir, prefix){
    
    # make result directory if not exist
    if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }

    # save seurat object
    saveRDS(seurat_obj, file=file.path(result_dir, paste0(prefix,"object",".rds")))
    # save metadata
    write.csv(seurat_obj[[]], file=file.path(result_dir, paste0(prefix,"metadata",".csv")), row.names=TRUE)
    # save stats
    stats <- paste0("cells: ",ncol(seurat_obj),"\nfeatures: ",nrow(seurat_obj))
    write(stats, file=file.path(result_dir, paste0(prefix,"stats",".txt")))
}


# extended ggsave
ggsave_new <- function(filename, results_path, plot, width=5, height=5){
    # make result directory if not exist
    if (!dir.exists(results_path)){
        dir.create(results_path, recursive = TRUE)
    }
    
    for (format in c('png')){
        ggsave(
          paste0(filename,'.',format),
          plot = plot,
          device = format,
          path = file.path(results_path),
          scale = 1,
          dpi = 300,
            width = width,
            height = height,
          limitsize = FALSE,
            units="in"
        )
    }
}