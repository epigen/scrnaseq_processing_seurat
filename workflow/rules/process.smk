# load & generate Seurat object per sample including metadata
rule prepare:
    input:
        get_sample_paths,
    output:
        sample_object = os.path.join(result_path,'batch_{sample}','prep_object.rds'),
        metadata = report(os.path.join(result_path,'batch_{sample}','prep_metadata.csv'), 
                          caption="../report/metadata_sample.rst", 
                          category="{}_scrnaseq_processing_seurat".format(config["project_name"]), 
                          subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","prepare_{sample}.log"),
    params:
        partition=config.get("partition"),
        sample = lambda w: "{}".format(w.sample),
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        crispr_umi_threshold = config["crispr_umi_threshold"],
        grna_regex = config["grna_regex"],
        percentage_regex = config["percentage_regex"],
        metadata_eval = config["metadata_eval"],
    script:
        "../scripts/prepare.R"

# merge into one dataset
rule merge:
    input:
        expand(os.path.join(result_path,'batch_{sample}','prep_object.rds'), sample=annot.index.tolist()),
        config["extra_metadata"] if config["extra_metadata"]!="" else [],
    output:
        merged_object = os.path.join(result_path,'merged','RAW_object.rds'),
        metadata = report(os.path.join(result_path,'merged','RAW_metadata.csv'), 
                          caption="../report/metadata_merged.rst", 
                          category="{}_scrnaseq_processing_seurat".format(config["project_name"]), 
                          subcategory="merged"),
    resources:
        mem_mb=config.get("mem", "16000"),
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
        merged_object = os.path.join(result_path,'merged','RAW_object.rds'),
    output:
        split_object = os.path.join(result_path,'{split}','RAW_object.rds'),
        metadata = report(os.path.join(result_path,'{split}','RAW_metadata.csv'), 
                          caption="../report/metadata.rst", 
                          category="{}_scrnaseq_processing_seurat".format(config["project_name"]), 
                          subcategory="{split}"),
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
        raw_object = os.path.join(result_path,'{split}','RAW_object.rds'),
    output:
        filtered_object = os.path.join(result_path,'{split}','FILTERED_object.rds'),
        metadata = report(os.path.join(result_path,'{split}','FILTERED_metadata.csv'), 
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
        step = "FILTERED",
        filter_expression = config["filter_expression"],
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/filter.R"
        
        
# save counts as CSV of seurat object
rule save_counts:
    input:
        seurat_object = os.path.join(result_path,'{split}','{step}_object.rds'),
    output:
        counts = os.path.join(result_path,'{split}','{step}_RNA.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
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