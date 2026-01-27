
#### load libraries & utility function 
library("Seurat")
library("sctransform")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@input[["utils_path"]])

# Set the future.globals.maxSize option to  90% of allocated memory, converting MB to bytes.
mem_mb <- as.numeric(snakemake@resources[["mem_mb"]])
options(future.globals.maxSize = (mem_mb * 0.9) * 1024^2)

# inputs
filtered_object_path <- snakemake@input[["filtered_object"]]
s_phase_genes_path <- snakemake@input[["s_phase_genes"]]
g2m_phase_genes_path <- snakemake@input[["g2m_phase_genes"]]


# outputs
norm_object_path <- snakemake@output[["norm_object"]]

# parameters
confounders <- snakemake@params[["confounders"]]
min_cells_per_gene <- snakemake@params[["min_cells_per_gene"]]
variable_features_n <- snakemake@params[["variable_features_n"]]
s_phase_genes <- snakemake@params[["s_phase_genes"]]
g2m_phase_genes <- snakemake@params[["g2m_phase_genes"]]
module_gene_lists <- snakemake@params[["module_gene_lists"]]
# 'flags' for modalities
ab_flag <- snakemake@params[['ab_flag']]
crispr_flag <- snakemake@params[['crispr_flag']]
custom_flag <- snakemake@params[['custom_flag']]


### load filtered data
filtered_object <- readRDS(file = file.path(filtered_object_path))

# load cell cycle scoring genes
if (s_phase_genes!=""){
    if (s_phase_genes=="tirosh2015"){
        s_genes <- cc.genes$s.genes
        g2m_genes <- cc.genes$g2m.genes
    }else{
        s_genes <- scan(file.path(s_phase_genes), character())
        g2m_genes <- scan(file.path(g2m_phase_genes), character())
    }
}

# load cell scoring gene list
gene_lists <- list()
for (gene_list_name in names(module_gene_lists)){
    tmp_genes <- scan(file.path(module_gene_lists[gene_list_name]), character())
    if(length(tmp_genes)>0){
        gene_lists[[gene_list_name]] <- tmp_genes
    }
}

### Normalization of all assays

# set parameters for HVF selection
if(variable_features_n==0){
    variable.features.n <- NULL
    variable.features.rv.th <- NULL
    return.only.var.genes <- FALSE
} else if(variable_features_n=="auto"){
    variable.features.n <- NULL
    variable.features.rv.th <- 1.3
    return.only.var.genes <- TRUE
} else {
    variable.features.n <- variable_features_n
    variable.features.rv.th <- NULL
    return.only.var.genes <- TRUE
}

# run SCTransform to normalize (& correct ie regress confounders) on RNA assay
norm_object <- SCTransform(filtered_object,
                           vars.to.regress = confounders,
                           assay = 'RNA',
                           new.assay.name = 'SCT',
                           verbose = TRUE,
                           method="glmGamPoi",
                           variable.features.n=variable.features.n,
                           variable.features.rv.th=variable.features.rv.th,
                           return.only.var.genes=return.only.var.genes,
                           min_cells = min_cells_per_gene,
                           seed.use = 42,
                           vst.flavor = "v2"
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


### Cell Scoring on normalized data (ie assay="SCT", slot="data")

# extract all features from the object
data_features <- rownames(GetAssayData(norm_object, slot = "data", assay = "SCT"))

# Cell Cycle scoring with Seurat function
# (presumably) running on SCT assay, as it is the default Assay post normalization
if (s_phase_genes!=""){
    norm_object <- CellCycleScoring(object = norm_object, s.features = s_genes, g2m.features = g2m_genes, search=TRUE) 
}

# Cell Scoring by gene list
for (gene_list_name in names(gene_lists)){
    
    # skip if no features in the data, to avoid Error
    if(length(intersect(gene_lists[[gene_list_name]], data_features))<3){
        print(paste0("Less than 3 of ",gene_list_name," features are present in the data."))
        next
    }

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

# higlhy variable genes (HVG)
if(length(confounders)==0){
    # get higlhy variable genes (HVG)
    HVG_df <- HVFInfo(object = norm_object[["SCT"]], selection.method = "sct")
    # save highly variable genes
    write(rownames(HVG_df)[order(-HVG_df$residual_variance)], file.path(dirname(file.path(norm_object_path)),"highly_variable_genes.txt"))
    fwrite(as.data.frame(HVG_df), file=file.path(dirname(file.path(norm_object_path)),"highly_variable_genes.csv"), row.names=TRUE)
}

# save all corrected data
save_seurat_object(seurat_obj=norm_object, result_dir=dirname(file.path(norm_object_path)))

