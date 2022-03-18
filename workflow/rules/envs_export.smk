rule env_export_seurat:
    output:
        os.path.join(config["result_path"],'envs','seurat.yaml'),
    conda:
        "../envs/seurat.yaml"
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_seurat.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        conda env export > {output}
        """