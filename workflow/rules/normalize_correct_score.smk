
# normalize data
rule normalize:
    input:
        filtered_object = os.path.join(result_path,'{split}','FILTERED_object.rds'),
    output:
        norm_object = os.path.join(result_path,'{split}','NORMALIZED_object.rds'),
        metadata = report(os.path.join(result_path,'{split}','NORMALIZED_metadata.csv'), 
                          caption="../report/metadata.rst", 
                          category="{}_scrnaseq_processing_seurat".format(config["project_name"]), 
                          subcategory="{split}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","filter_{split}.log"),
    params:
        partition=config.get("partition"),
        step = "NORMALIZED",
        min_cells_per_gene = config["min_cells_per_gene"],
        confounders = [],
        module_gene_lists = config["module_gene_lists"],
        cell_cycle = config["cell_cycle"],
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/sctransform_cellScore.R"
        
        
        
# correct data
if len(config["variables_to_regress"])>0:
    rule correct:
        input:
            filtered_object = os.path.join(result_path,'{split}','NORMALIZED_object.rds'),
        output:
            norm_object = os.path.join(result_path,'{split}','CORRECTED_object.rds'),
            metadata = report(os.path.join(result_path,'{split}','CORRECTED_metadata.csv'), 
                              caption="../report/metadata.rst", 
                              category="{}_scrnaseq_processing_seurat".format(config["project_name"]), 
                              subcategory="{split}"),
        resources:
            mem_mb=config.get("mem", "16000"),
        threads: config.get("threads", 1)
        conda:
            "../envs/seurat.yaml"
        log:
            os.path.join("logs","rules","filter_{split}.log"),
        params:
            partition=config.get("partition"),
            step = "CORRECTED",
            min_cells_per_gene = config["min_cells_per_gene"],
            confounders = config["variables_to_regress"],
            module_gene_lists = config["module_gene_lists"],
            cell_cycle = config["cell_cycle"],
            ab_flag = config["modality_flags"]['Antibody_Capture'],
            crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
            custom_flag = config["modality_flags"]['Custom'],
        script:
            "../scripts/sctransform_cellScore.R"