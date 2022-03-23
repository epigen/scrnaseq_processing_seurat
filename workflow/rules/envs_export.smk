# one rule per used conda environment to document the exact versions and builds of the used software

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
        
        
rule env_export_inspectdf:
    output:
        os.path.join(config["result_path"],'envs','inspectdf.yaml'),
    conda:
        "../envs/inspectdf.yaml"
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_inspectdf.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        conda env export > {output}
        """