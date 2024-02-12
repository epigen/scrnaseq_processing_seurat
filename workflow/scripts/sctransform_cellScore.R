
#### load libraries & utility function 
library("Seurat")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
filtered_object_path <- snakemake@input[["filtered_object"]]

# outputs
norm_object_path <- snakemake@output[["norm_object"]]

# parameters
confounders <- snakemake@params[["confounders"]]
min_cells_per_gene <- snakemake@params[["min_cells_per_gene"]]

module_gene_lists <- snakemake@params[["module_gene_lists"]]
cell_cycle <- snakemake@params[["cell_cycle"]]

ab_flag <- snakemake@params[["ab_flag"]]
crispr_flag <- snakemake@params[["crispr_flag"]]
custom_flag <- snakemake@params[["custom_flag"]]


### load filtered data
filtered_object <- readRDS(file = file.path(filtered_object_path))

# load cell cycle scoring genes
if (cell_cycle['s_phase_genes']=="tirosh2015"){
    s_genes <- cc.genes$s.genes
    g2m_genes <- cc.genes$g2m.genes
}else{
    s_genes <- scan(file.path(cell_cycle['s_phase_genes']), character())
    g2m_genes <- scan(file.path(cell_cycle['g2m_phase_genes']), character())
}

# load cell scoring gene list
gene_lists <- list()
for (gene_list_name in names(module_gene_lists)){
    gene_lists[[gene_list_name]] <- scan(file.path(module_gene_lists[gene_list_name]), character())
}


### Normalization of all assays

# run SCTransform to normalize (& correct ie regress confounders) on RNA assay
norm_object <- SCTransform(filtered_object,
                           vars.to.regress = confounders,
                           assay = 'RNA',
                           new.assay.name = 'SCT',
                           verbose = TRUE,
                           method="glmGamPoi",
                           variable.features.n=NULL,
                           return.only.var.genes=FALSE,
                           min_cells = min_cells_per_gene,
                           seed.use = 42
                          )

# normalize AB data (margin=2 -> normalization across cells)
if(ab_flag!=''){
    norm_object <- NormalizeData(norm_object, normalization.method = "CLR", margin = 2, assay = ab_flag)
    }

# normalize CRISPR data (margin=2 -> normalization across cells)
if(crispr_flag!=''){
    norm_object <- NormalizeData(norm_object, normalization.method = "CLR", margin = 2, assay = crispr_flag)
}

# normalize Custom data (margin=2 -> normalization across cells)
if(custom_flag!=''){
    norm_object <- NormalizeData(norm_object, normalization.method = "CLR", margin = 2, assay = custom_flag)
}


# get higlhy variable genes (HVG)
HVG_df <- HVFInfo(object = norm_object, assay = "SCT")


### Cell Scoring on normalized data (ie assay="SCT", slot="data")

# Cell Cycle scoring with integrated Seurat function
# (presumably) running on SCT assay, as it is the default Assay post normalization
if (cell_cycle['s_phase_genes']!=""){
    norm_object <- CellCycleScoring(object = norm_object, s.features = s_genes, g2m.features = g2m_genes, search=TRUE) 
}

# Cell Scoring by gene list
for (gene_list_name in names(gene_lists)){
    norm_object <- AddModuleScore(
        object=norm_object,
        features=list(gene_lists[[gene_list_name]]),
        pool = NULL,
        nbin = 24,
        ctrl = 100,
        k = FALSE,
        assay = "SCT",
        name = gene_list_name,
        seed = 42,
        search = TRUE, #Search for symbol synonyms for features in features that don't match features in object (ie gene symbol synonyms) using HUGO Gene Nomenclature Committee (HGNC)
    )
    
    # remove post-fix '1' from metadata column
    norm_object[[gene_list_name]] <- norm_object[[paste0(gene_list_name,'1')]]
    norm_object[[paste0(gene_list_name,'1')]] <- NULL
}


### save data
    
# save highly variable genes
write(rownames(HVG_df), file.path(dirname(file.path(norm_object_path)),"highly_variable_genes.txt"))
#     write.csv(HVG_df, file=file.path(dirname(file.path(norm_object_path)),"highly_variable_genes.csv"), row.names=TRUE)
fwrite(as.data.frame(HVG_df), file=file.path(dirname(file.path(norm_object_path)),"highly_variable_genes.csv"), row.names=TRUE)

# save all corrected data
save_seurat_object(seurat_obj=norm_object, result_dir=dirname(file.path(norm_object_path)))

