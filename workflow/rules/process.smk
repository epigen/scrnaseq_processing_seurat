
# load & generate Seurat object per sample including metadata
rule prepare:
    input:
        sample_dir = lambda w: annot.loc["{}".format(w.sample),'path'], # get_sample_path,
        metadata = lambda w: [] if pd.isna(annot.loc["{}".format(w.sample),'metadata']) else annot.loc["{}".format(w.sample),'metadata'],
        utils_path = workflow.source_path("../scripts/utils.R"),
    output:
        sample_object = os.path.join(result_path,'batch__{sample}','PREP','object.rds'),
        metadata = report(os.path.join(result_path,'batch__{sample}','PREP','metadata.csv'), 
                          caption="../report/metadata_sample.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{sample}",
                          labels={
                              "step": "PREP",
                              "type": "CSV",
                              "category": "All",
                              "list": "Metadata",
                              "feature": "",
                                }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","PREP_{sample}.log"),
    params:
        # metadata = lambda w: "" if pd.isna(annot.loc["{}".format(w.sample),'metadata']) else annot.loc["{}".format(w.sample),'metadata'],
    script:
        "../scripts/prepare.R"

# merge into one dataset
rule merge:
    input:
        samples = expand(os.path.join(result_path,'batch__{sample}','PREP','object.rds'), sample=annot.index.tolist()),
        extra_metadata = config["extra_metadata"] if config["extra_metadata"]!="" else [],
        utils_path = workflow.source_path("../scripts/utils.R"),
    output:
        merged_object = os.path.join(result_path,'merged','PREP','object.rds'),
        metadata = report(os.path.join(result_path,'merged','PREP','metadata.csv'), 
                          caption="../report/metadata_merged.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="merged",
                          labels={
                              "step": "RAW",
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
        os.path.join("logs","rules","merge.log"),
    script:
        "../scripts/merge.R"

# split into subsets by metadata
rule split:
    wildcard_constraints:
        split=r"(?!merged$|batch__)[^/]+"
    input:
        merged_object = os.path.join(result_path,'merged','PREP','object.rds'),
        utils_path = workflow.source_path("../scripts/utils.R"),
    output:
        split_object = os.path.join(result_path,'{split}','RAW','object.rds'),
        metadata = report(os.path.join(result_path,'{split}','RAW','metadata.csv'), 
                          caption="../report/metadata.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "step": "RAW",
                              "type": "CSV",
                              "category": "All",
                              "list": "Metadata",
                              "feature": "",
                                }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","split_{split}.log"),
    script:
        "../scripts/split.R"

# filter cells according to metadata
rule filter_cells:
    input:
        raw_object = os.path.join(result_path,'{split}','RAW','object.rds'),
        utils_path = workflow.source_path("../scripts/utils.R"),
    output:
        filtered_object = os.path.join(result_path,'{split}','FILTERED','object.rds'),
        metadata = report(os.path.join(result_path,'{split}','FILTERED','metadata.csv'), 
                          caption="../report/metadata.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "step": "FILTERED",
                              "type": "CSV",
                              "category": "All",
                              "list": "Metadata",
                              "feature": "",
                                }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","filter_{split}.log"),
    params:
        step = "FILTERED",
        filter_expression = config["filter_expression"],
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/filter.R"

# pseudobulk cells into samples
rule pseudobulk:
    input:
        filtered_object = os.path.join(result_path,'{split}','FILTERED','object.rds'),
        utils_path = workflow.source_path("../scripts/utils.R"),
    output:
        pseudobulk_counts = os.path.join(result_path,'{split}','PSEUDOBULK','RNA.csv'),
        metadata = report(os.path.join(result_path,'{split}','PSEUDOBULK','metadata.csv'), 
                          caption="../report/metadata.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "step": "PSEUDOBULK",
                              "type": "CSV",
                              "category": "All",
                              "list": "Metadata",
                              "feature": "",
                                }),
        cell_count_plot = report(os.path.join(result_path,'{split}','PSEUDOBULK','cell_count_histogram_density.png'), 
                          caption="../report/pseudobulk_cell_count.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{split}",
                          labels={
                              "step": "PSEUDOBULK",
                              "type": "Histogram",
                              "category": "All",
                              "list": "Cell count",
                              "feature": "",
                                }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","pseudobulk_{split}.log"),
    script:
        "../scripts/pseudobulk.R"

# save counts as CSV of seurat object
rule save_counts:
    input:
        seurat_object = os.path.join(result_path,'{split}','{step}','object.rds'),
    output:
        rna = os.path.join(result_path,'{split}','{step}','RNA.csv'),
        ab = [] if (config["modality_flags"]['Antibody_Capture'] == "") else os.path.join(result_path, '{split}', '{step}', f'{config["modality_flags"]["Antibody_Capture"]}.csv'),
        crispr = [] if (config["modality_flags"]['CRISPR_Guide_Capture'] == "") else os.path.join(result_path, '{split}', '{step}', f'{config["modality_flags"]["CRISPR_Guide_Capture"]}.csv'),
        custom = [] if (config["modality_flags"]['Custom'] == "") else os.path.join(result_path, '{split}', '{step}', f'{config["modality_flags"]["Custom"]}.csv'),
    resources:
        mem_mb = lambda wc, attempt: attempt*int(config.get("mem", "16000")),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","save_counts_{split}_{step}.log"),
    params:
        step = lambda w: "{}".format(w.step),
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/save_counts.R"