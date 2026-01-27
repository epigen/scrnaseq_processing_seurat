
# visualize metadata
rule metadata_plots:
    input:
        metadata = os.path.join(result_path,'{split}','{step}','metadata.csv'),
        utils_path = workflow.source_path("../scripts/utils.R"),
    output:
        metadata_plots = report(directory(os.path.join(result_path,'{split}','{step}','plots','metadata')),
                                patterns=["{datatype}.png"],
                                caption="../report/metadata_vis.rst", 
                                category="{}_{}".format(config["project_name"], module_name), 
                                subcategory="{split}",
                                labels={
                                "step": "{step}",
                                "type": "{datatype}",
                                "category": "",
                                "list": "Metadata",
                                "feature": "",
                                }),
        metadata_stats = directory(os.path.join(result_path,'{split}','{step}','stats')),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/inspectdf.yaml"
    log:
        os.path.join("logs","rules","metadata_plot_{split}_{step}.log"),
    script:
        "../scripts/metadata_plot.R"

# make Seurat plots (ridge, violin, dot plot, heatmap)
rule seurat_plots:
    input:
        norm_object = os.path.join(result_path,'{split}','{step}','object.rds'),
        vis_gene_lists = get_vis_gene_lists, # for tracking inputs
        utils_path = workflow.source_path("../scripts/utils.R"),
    output:
        plot_dir = report(
            directory(os.path.join(result_path,'{split}','{step}','plots','{plot_type}','{category}','{feature_list}')),
            patterns=["{feature}.png"],
            caption="../report/seurat_plot.rst",
            category="{}_{}".format(config["project_name"], module_name),
            subcategory="{split}",
            labels={
                "step": "{step}",
                "type": "{plot_type}",
                "category": "{category}",
                "list": "{feature_list}",
                "feature": "{feature}",
            }),
    params:
        vis_gene_lists = config["vis_gene_lists"],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","plot_{split}_{step}_{plot_type}_{category}_{feature_list}.log"),
    script:
        "../scripts/seurat_plots.R"
