
# normalize data
rule normalize:
    input:
        filtered_object = os.path.join(config["result_path"],'{split}','counts','FILTERED_object.rds'),
    output:
        norm_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
        metadata = report(os.path.join(config["result_path"],'{split}','counts','NORMALIZED_metadata.csv'), caption="../report/metadata.rst", category="processing_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","filter_{split}.log"),
    params:
        partition=config.get("partition"),
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
            filtered_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
        output:
            norm_object = os.path.join(config["result_path"],'{split}','counts','CORRECTED_object.rds'),
            metadata = report(os.path.join(config["result_path"],'{split}','counts','CORRECTED_metadata.csv'), caption="../report/metadata.rst", category="processing_{}".format(config["project_name"]), subcategory="{split}"),
        resources:
            mem=config.get("mem", "16G"),
        threads: config.get("threads", 1)
        conda:
            "../envs/seurat.yaml"
        log:
            os.path.join("logs","rules","filter_{split}.log"),
        params:
            partition=config.get("partition"),
            min_cells_per_gene = config["min_cells_per_gene"],
            confounders = config["variables_to_regress"],
            module_gene_lists = config["module_gene_lists"],
            cell_cycle = config["cell_cycle"],
            ab_flag = config["modality_flags"]['Antibody_Capture'],
            crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
            custom_flag = config["modality_flags"]['Custom'],
        script:
            "../scripts/sctransform_cellScore.R"