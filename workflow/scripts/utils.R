### Utility Functions

# save Seurat object results
save_seurat_object <- function (seurat_obj, result_dir, prefix, rna_flag="RNA", ab_flag, crispr_flag, custom_flag, slot="counts", saveCounts){
    
    # make result directory if not exist
    if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }
    
    # to avoid multiple saving of the seurat object & metadata post normalization and correction
    if (!grepl('counts_', prefix, fixed = TRUE)){
        # save seurat object
        saveRDS(seurat_obj, file=file.path(result_dir, paste0(prefix,"object",".rds")))
        # save metadata
        write.csv(seurat_obj[[]], file=file.path(result_dir, paste0(prefix,"metadata",".csv")), row.names=TRUE)
    }

    if(saveCounts==1){
        # save slot data of RNA assay
        write.csv(GetAssayData(object = seurat_obj, slot = slot, assay = rna_flag), file=file.path(result_dir, paste0(prefix,rna_flag,'.csv')), row.names=TRUE)

        # normalized values of all other modalities are in the "data" slot
        if (slot=="scale.data"){
            slot <- "data"
        }

        # save slot data of all assays
        if(ab_flag!=''){
            write.csv(GetAssayData(object = seurat_obj, slot = slot, assay = ab_flag), file=file.path(result_dir, paste0(prefix,ab_flag,".csv")), row.names=TRUE)
        }

        if(crispr_flag!=''){
            write.csv(GetAssayData(object = seurat_obj, slot = slot, assay = crispr_flag), file=file.path(result_dir, paste0(prefix,crispr_flag,'.csv')), row.names=TRUE)
        }

        if(custom_flag!=''){
            write.csv(GetAssayData(object = seurat_obj, slot = slot, assay = custom_flag), file=file.path(result_dir, paste0(prefix,custom_flag,'.csv')), row.names=TRUE)
        }
    }
    
}