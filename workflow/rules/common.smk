##### utility functions #####

def get_sample_paths(wildcards):
    return [annot.loc[wildcards.sample,'path'], annot.loc[wildcards.sample,'metadata']]



# save counts as CSV of seurat object
rule save_counts:
    input:
        seurat_object = os.path.join(config["result_path"],'{split}','counts','{step}_object.rds'),
    output:
        counts = os.path.join(config["result_path"],'{split}','counts','{step}_RNA.csv'),
    resources:
        mem=config.get("mem", "16G"),
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