
# normalize data
rule normalize:
    input:
        filtered_object = os.path.join(config["result_path"],'{split}','counts','FILTERED_object.rds'),
    output:
        norm_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","filter_{split}.log"),
    params:
        partition=config.get("partition"),
        saveCounts = config["saveCounts"],
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
            filtered_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
        output:
            norm_object = os.path.join(config["result_path"],'{split}','counts','CORRECTED_object.rds'),
        resources:
            mem=config.get("mem", "16G"),
        threads: config.get("threads", 1)
        conda:
            "../envs/seurat.yaml"
        log:
            os.path.join("logs","rules","filter_{split}.log"),
        params:
            partition=config.get("partition"),
            saveCounts = config["saveCounts"],
            confounders = config["variables_to_regress"],
            module_gene_lists = config["module_gene_lists"],
            cell_cycle = config["cell_cycle"],
            ab_flag = config["modality_flags"]['Antibody_Capture'],
            crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
            custom_flag = config["modality_flags"]['Custom'],
        script:
            "../scripts/sctransform_cellScore.R"