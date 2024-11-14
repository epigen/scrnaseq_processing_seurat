
# normalize data
rule normalize:
    input:
        filtered_object = os.path.join(result_path,'{split}','FILTERED','object.rds'),
    output:
        norm_object = os.path.join(result_path,'{split}','NORMALIZED','object.rds'),
        metadata = report(os.path.join(result_path,'{split}','NORMALIZED','metadata.csv'), 
                          caption="../report/metadata.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "step": "NORMALIZED",
                              "type": "CSV",
                              "category": "All",
                              "list": "Metadata",
                              "feature": "",
                                }),
    resources:
        mem_mb = lambda wc, attempt: attempt*int(config.get("mem", "16000")),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","filter_{split}.log"),
    params:
        min_cells_per_gene = config["min_cells_per_gene"],
        variable_features_n = config["variable_features_n"],
        confounders = [],
        module_gene_lists = config["module_gene_lists"],
        cell_cycle = config["cell_cycle"],
    script:
        "../scripts/sctransform_cellScore.R"
        
        
        
# correct data
rule correct:
    input:
        # NORMALIZED object as input as only it contains post-normalization calculated scores to be regressed out e.g., cell-cycle scores
        filtered_object = os.path.join(result_path,'{split}','NORMALIZED','object.rds'), 
    output:
        norm_object = os.path.join(result_path,'{split}','CORRECTED','object.rds'),
        metadata = report(os.path.join(result_path,'{split}','CORRECTED','metadata.csv'), 
                          caption="../report/metadata.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "step": "CORRECTED",
                              "type": "CSV",
                              "category": "All",
                              "list": "Metadata",
                              "feature": "",
                          }),
    resources:
        mem_mb = lambda wc, attempt: attempt*int(config.get("mem", "16000")),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","correct_{split}.log"),
    params:
        min_cells_per_gene = config["min_cells_per_gene"],
        variable_features_n = config["variable_features_n"],
        confounders = config["variables_to_regress"],
        module_gene_lists = config["module_gene_lists"],
        cell_cycle = config["cell_cycle"],
    script:
        "../scripts/sctransform_cellScore.R"