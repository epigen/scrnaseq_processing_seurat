
#### load libraries & utility function 
library("Seurat")
library("data.table")

# inputs
object_path <- snakemake@input[["seurat_object"]]

# outputs
counts_path <- snakemake@output[["counts"]]

# parameters
step <- snakemake@params[["step"]]
# 'flags' for modalities
ab_flag <- snakemake@params[["ab_flag"]]
crispr_flag <- snakemake@params[["crispr_flag"]]
custom_flag <- snakemake@params[["custom_flag"]]

result_dir <- file.path(dirname(counts_path))

### load seurat data
seurat_object <- readRDS(file = file.path(object_path))

### save count data

# save slot data of RNA assay depending on processing step
if (step=="RAW" | step=="FILTERED"){
    fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = "counts", assay = "RNA")), file=file.path(result_dir,'RNA.csv'), row.names=TRUE)
} else if (step=="NORMALIZED"){
    fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = "counts", assay = "SCT")), file=file.path(result_dir,'RNA_counts.csv'), row.names=TRUE)
    fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = "data", assay = "SCT")), file=file.path(result_dir,'RNA.csv'), row.names=TRUE)
    fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = "scale.data", assay = "SCT")), file=file.path(result_dir,'RNA_scaled.csv'), row.names=TRUE)
} else if (step=="CORRECTED") {
    fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = "scale.data", assay = "SCT")), file=file.path(result_dir,'RNA.csv'), row.names=TRUE)
}

# save slot data of all other assays, if step=="CORRECTED" then skipped (as other assays are not corrected)
if (step!="CORRECTED"){
    # normalized values of all other modalities are in the "data" slot
    if (step=="NORMALIZED"){
        slot <- "data"
    }else{
        slot <- "counts"
    }

    if(ab_flag!=''){
        fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = slot, assay = ab_flag)), file=file.path(result_dir, paste0(ab_flag,".csv")), row.names=TRUE)
    }

    if(crispr_flag!=''){
        fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = slot, assay = crispr_flag)), file=file.path(result_dir, paste0(crispr_flag,'.csv')), row.names=TRUE)
    }

    if(custom_flag!=''){
        fwrite(as.data.frame(GetAssayData(object = seurat_object, slot = slot, assay = custom_flag)), file=file.path(result_dir, paste0(custom_flag,'.csv')), row.names=TRUE)
    }
}
