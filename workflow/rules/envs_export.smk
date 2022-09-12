# one rule per used conda environment to document the exact versions and builds of the used software        
rule env_export:
    output:
        report(os.path.join(config["result_path"],'envs','scrnaseq_processing_seurat','{env}.yaml'),
                      caption="../report/software.rst", 
                      category="Software", 
                      subcategory="{}_scrnaseq_processing_seurat".format(config["project_name"])
                     ),
    conda:
        "../envs/{env}.yaml"
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_{env}.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        conda env export > {output}
        """
        
# add configuration files to report        
rule config_export:
    output:
        configs = report(os.path.join(config["result_path"],'configs','scrnaseq_processing_seurat','{}_config.yaml'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_scrnaseq_processing_seurat".format(config["project_name"])
                        )
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","config_export.log"),
    params:
        partition=config.get("partition"),
    run:
        with open(output["configs"], 'w') as outfile:
            yaml.dump(config, outfile)

# export used annotation file for documentation and reproducibility            
rule annot_export:
    output:
        configs = report(config["sample_annotation"], 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_scrnaseq_processing_seurat".format(config["project_name"])
                        )
    resources:
        mem_mb=config.get("mem_small", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","annot_export.log"),
    params:
        partition=config.get("partition"),

# # export used gene lists for documentation and reproducibility -> problems: 1) if error occurs they get deleted as they are an output. 2) if mutliple modules use the same gene lists in a project, the report can not be generated due to an AmbiguousRuleException.
# rule gene_list_export:
#     output:
#         gene_lists = report(set(list(config["module_gene_lists"].values())+list(config["vis_gene_lists"].values())), 
#                             caption="../report/gene_lists.rst", 
#                             category="Configuration", 
#                             subcategory="{}_scrnaseq_processing_seurat".format(config["project_name"])
#                            ),
#     resources:
#         mem_mb=config.get("mem", "16000"),
#     threads: config.get("threads", 1)
#     log:
#         os.path.join("logs","rules","gene_list_export.log"),
#     params:
#         partition=config.get("partition"),

        

        