# load & generate Seurat object per sample including metadata
rule prepare:
    input:
        get_sample_paths,
    output:
        sample_object = os.path.join(config["result_path"],'{sample}','counts','RAW_object.rds'),
        metadata = report(os.path.join(config["result_path"],'{sample}','counts','RAW_metadata.csv'), caption="../report/metadata_sample.rst", category="processing_{}".format(config["project_name"]), subcategory="{sample}"),
    resources:
        mem=config.get("mem", "16G"),
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
        expand(os.path.join(config["result_path"],'{sample}','counts','RAW_object.rds'), sample=annot.index.tolist()),
    output:
        merged_object = os.path.join(config["result_path"],'merged','counts','RAW_object.rds'),
        metadata = report(os.path.join(config["result_path"],'merged','counts','RAW_metadata.csv'), caption="../report/metadata_merged.rst", category="processing_{}".format(config["project_name"]), subcategory="merged"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","merge.log"),
    params:
        partition=config.get("partition"),
        project_name = config["project_name"],
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/merge.R"

# split into subsets by metadata
rule split:
    input:
        merged_object = os.path.join(config["result_path"],'merged','counts','RAW_object.rds'),
    output:
        split_object = os.path.join(config["result_path"],'{split}','counts','RAW_object.rds'),
        metadata = report(os.path.join(config["result_path"],'{split}','counts','RAW_metadata.csv'), caption="../report/metadata.rst", category="processing_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","split_{split}.log"),
    params:
        partition=config.get("partition"),
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
        raw_object = os.path.join(config["result_path"],'{split}','counts','RAW_object.rds'),
    output:
        filtered_object = os.path.join(config["result_path"],'{split}','counts','FILTERED_object.rds'),
        metadata = report(os.path.join(config["result_path"],'{split}','counts','FILTERED_metadata.csv'), caption="../report/metadata.rst", category="processing_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","filter_{split}.log"),
    params:
        partition=config.get("partition"),
        filter_expression = config["filter_expression"],
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
    script:
        "../scripts/filter.R"