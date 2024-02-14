# load & generate Seurat object per sample including metadata
rule prepare:
    input:
        get_sample_paths,
    output:
        sample_object = os.path.join(result_path,'batch_{sample}','PREP','object.rds'),
        metadata = report(os.path.join(result_path,'batch_{sample}','PREP','metadata.csv'), 
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
        partition=config.get("partition"),
        metadata = lambda w: "" if pd.isna(annot.loc["{}".format(w.sample),'metadata']) else annot.loc["{}".format(w.sample),'metadata'],
    script:
        "../scripts/prepare.R"

# merge into one dataset
rule merge:
    input:
        expand(os.path.join(result_path,'batch_{sample}','PREP','object.rds'), sample=annot.index.tolist()),
        config["extra_metadata"] if config["extra_metadata"]!="" else [],
    output:
        merged_object = os.path.join(result_path,'merged','RAW','object.rds'),
        metadata = report(os.path.join(result_path,'merged','RAW','metadata.csv'), 
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
#         mem_mb=config.get("mem", "16000"),
        mem_mb = lambda wc, attempt: attempt*int(config.get("mem", "16000")),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","merge.log"),
    params:
        partition=config.get("partition"),
        step = "RAW",
        project_name = config["project_name"],
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        extra_metadata = config["extra_metadata"],
    script:
        "../scripts/merge.R"

# split into subsets by metadata
rule split:
    input:
        merged_object = os.path.join(result_path,'merged','RAW','object.rds'),
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
    params:
        partition=config.get("partition"),
        step = "RAW",
        result_dir = lambda w, input: os.path.splitext(input[0])[0],
        split = lambda w: "{}".format(w.split),
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/split.R"

# filter cells according to metadata
rule filter_cells:
    input:
        raw_object = os.path.join(result_path,'{split}','RAW','object.rds'),
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
        partition=config.get("partition"),
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
    params:
        partition=config.get("partition"),
    script:
        "../scripts/pseudobulk.R"

# save counts as CSV of seurat object
rule save_counts:
    input:
        seurat_object = os.path.join(result_path,'{split}','{step}','object.rds'),
    output:
        counts = os.path.join(result_path,'{split}','{step}','RNA.csv'),
    resources:
#         mem_mb=config.get("mem", "16000"),
        mem_mb = lambda wc, attempt: attempt*int(config.get("mem", "16000")),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","save_counts_{split}_{step}.log"),
    params:
        partition=config.get("partition"),
        step = lambda w: "{}".format(w.step),
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/save_counts.R"