# visualize metadata
rule metadata_plots:
    input:
        metadata = os.path.join(config["result_path"],'{split}','counts','{step}_metadata.csv'),
    output:
        metadata_plots = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','{{step}}_metadata_{datatype}.png'),
                             datatype=['numerical','categorical','types']
                            ), caption="../report/metadata_vis.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/inspectdf.yaml"
    log:
        os.path.join("logs","rules","metadata_plot_{split}_{step}.log"),
    params:
        partition=config.get("partition"),
        step = lambda w: "{}".format(w.step),
    script:
        "../scripts/metadata_plot.R"

# make ridge plots
rule ridge_plot_normalized:
    input:
        norm_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
    output:
        ridge_plots_expression = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_ridge_plot_{category}_{gene_list}.png'),
                             category=config["vis_categories"],
                             gene_list=list(config["vis_gene_lists"].keys())
                            ), caption="../report/ridge_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        ridge_plots_Antibody_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_ridge_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Antibody_Capture"]) if (config["modality_flags"]["Antibody_Capture"]!="")and(len(config["vis_features"]["Antibody_Capture"])>0) else None, caption="../report/ridge_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        ridge_plots_CRISPR_Guide_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_ridge_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["CRISPR_Guide_Capture"]) if (config["modality_flags"]["CRISPR_Guide_Capture"]!="")and(len(config["vis_features"]["CRISPR_Guide_Capture"])>0) else None, caption="../report/ridge_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        ridge_plots_Custom = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_ridge_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Custom"]) if (config["modality_flags"]["Custom"]!="")and(len(config["vis_features"]["Custom"])>0) else None, caption="../report/ridge_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        ridge_plots_Metadata = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_ridge_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality="Metadata") if (len(config["vis_features"]["Metadata"])>0) else None, caption="../report/ridge_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","ridge_plot_{split}_NORMALIZED.log"),
    params:
        partition=config.get("partition"),
        step = "NORMALIZED",
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        vis_categories = config["vis_categories"],
        vis_gene_lists = config["vis_gene_lists"],
        ab_features = config["vis_features"]['Antibody_Capture'],
        crispr_features = config["vis_features"]['CRISPR_Guide_Capture'],
        custom_features = config["vis_features"]['Custom'],
        metadata_features = config["vis_features"]['Metadata'],
    script:
        "../scripts/ridge_plot.R"
        
rule ridge_plot_corrected:
    input:
        norm_object = os.path.join(config["result_path"],'{split}','counts','CORRECTED_object.rds'),
    output:
        ridge_plots_expression = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','CORRECTED_ridge_plot_{category}_{gene_list}.png'),
                             category=config["vis_categories"],
                             gene_list=list(config["vis_gene_lists"].keys())
                            ), caption="../report/ridge_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","ridge_plot_{split}_CORRECTED.log"),
    params:
        partition=config.get("partition"),
        step = "CORRECTED",
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        vis_categories = config["vis_categories"],
        vis_gene_lists = config["vis_gene_lists"],
        ab_features = config["vis_features"]['Antibody_Capture'],
        crispr_features = config["vis_features"]['CRISPR_Guide_Capture'],
        custom_features = config["vis_features"]['Custom'],
    script:
        "../scripts/ridge_plot.R"
        
        
# make violin plots
rule violin_plot_normalized:
    input:
        norm_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
    output:
        violin_plots_expression = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_violin_plot_{category}_{gene_list}.png'),
                             category=config["vis_categories"],
                             gene_list=list(config["vis_gene_lists"].keys())
                            ), caption="../report/violin_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        violin_plots_Antibody_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_violin_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Antibody_Capture"]) if (config["modality_flags"]["Antibody_Capture"]!="")and(len(config["vis_features"]["Antibody_Capture"])>0) else None, caption="../report/violin_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        violin_plots_CRISPR_Guide_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_violin_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["CRISPR_Guide_Capture"]) if (config["modality_flags"]["CRISPR_Guide_Capture"]!="")and(len(config["vis_features"]["CRISPR_Guide_Capture"])>0) else None, caption="../report/violin_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        violin_plots_Custom = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_violin_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Custom"]) if (config["modality_flags"]["Custom"]!="")and(len(config["vis_features"]["Custom"])>0) else None, caption="../report/violin_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        violin_plots_Metadata = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_violin_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality="Metadata") if (len(config["vis_features"]["Metadata"])>0) else None, caption="../report/violin_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","violin_plot_{split}_NORMALIZED.log"),
    params:
        partition=config.get("partition"),
        step = "NORMALIZED",
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        vis_categories = config["vis_categories"],
        vis_gene_lists = config["vis_gene_lists"],
        ab_features = config["vis_features"]['Antibody_Capture'],
        crispr_features = config["vis_features"]['CRISPR_Guide_Capture'],
        custom_features = config["vis_features"]['Custom'],
        metadata_features = config["vis_features"]['Metadata'],
    script:
        "../scripts/violin_plot.R"
        
rule violin_plot_corrected:
    input:
        norm_object = os.path.join(config["result_path"],'{split}','counts','CORRECTED_object.rds'),
    output:
        violin_plots_expression = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','CORRECTED_violin_plot_{category}_{gene_list}.png'),
                             category=config["vis_categories"],
                             gene_list=list(config["vis_gene_lists"].keys())
                            ), caption="../report/violin_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","violin_plot_{split}_CORRECTED.log"),
    params:
        partition=config.get("partition"),
        step = "CORRECTED",
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        vis_categories = config["vis_categories"],
        vis_gene_lists = config["vis_gene_lists"],
        ab_features = config["vis_features"]['Antibody_Capture'],
        crispr_features = config["vis_features"]['CRISPR_Guide_Capture'],
        custom_features = config["vis_features"]['Custom'],
    script:
        "../scripts/violin_plot.R"
        
# make heatmaps
rule heatmap_normalized:
    input:
        norm_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
    output:
        heatmaps_expression = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_heatmap_{category}_{gene_list}.png'),
                             category=config["vis_categories"],
                             gene_list=list(config["vis_gene_lists"].keys())+['HVG100']
                            ), caption="../report/heatmap.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        heatmaps_Antibody_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_heatmap_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Antibody_Capture"]) if (config["modality_flags"]["Antibody_Capture"]!="")and(len(config["vis_features"]["Antibody_Capture"])>0) else None, caption="../report/heatmap.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        heatmaps_CRISPR_Guide_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_heatmap_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["CRISPR_Guide_Capture"]) if (config["modality_flags"]["CRISPR_Guide_Capture"]!="")and(len(config["vis_features"]["CRISPR_Guide_Capture"])>0) else None, caption="../report/heatmap.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        heatmaps_Custom = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_heatmap_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Custom"]) if (config["modality_flags"]["Custom"]!="")and(len(config["vis_features"]["Custom"])>0) else None, caption="../report/heatmap.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","heatmap_{split}_NORMALIZED.log"),
    params:
        partition=config.get("partition"),
        step = "NORMALIZED",
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        vis_categories = config["vis_categories"],
        vis_gene_lists = config["vis_gene_lists"],
        ab_features = config["vis_features"]['Antibody_Capture'],
        crispr_features = config["vis_features"]['CRISPR_Guide_Capture'],
        custom_features = config["vis_features"]['Custom'],
    script:
        "../scripts/heatmap_plot.R"
        
rule heatmap_corrected:
    input:
        norm_object = os.path.join(config["result_path"],'{split}','counts','CORRECTED_object.rds'),
    output:
        heatmaps_expression = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','CORRECTED_heatmap_{category}_{gene_list}.png'),
                             category=config["vis_categories"],
                             gene_list=list(config["vis_gene_lists"].keys())+['HVG100']
                            ), caption="../report/heatmap.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","heatmap_{split}_CORRECTED.log"),
    params:
        partition=config.get("partition"),
        step = "CORRECTED",
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        vis_categories = config["vis_categories"],
        vis_gene_lists = config["vis_gene_lists"],
        ab_features = config["vis_features"]['Antibody_Capture'],
        crispr_features = config["vis_features"]['CRISPR_Guide_Capture'],
        custom_features = config["vis_features"]['Custom'],
    script:
        "../scripts/heatmap_plot.R"        
        
# make dot plots
rule dot_plot:
    input:
        norm_object = os.path.join(config["result_path"],'{split}','counts','NORMALIZED_object.rds'),
    output:
        dot_plots_expression = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_dot_plot_{category}_{gene_list}.png'),
                             category=config["vis_categories"],
                             gene_list=list(config["vis_gene_lists"].keys())
                            ), caption="../report/dot_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        dot_plots_Antibody_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_dot_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Antibody_Capture"]) if (config["modality_flags"]["Antibody_Capture"]!="")and(len(config["vis_features"]["Antibody_Capture"])>0) else None, caption="../report/dot_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        dot_plots_CRISPR_Guide_Capture = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_dot_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["CRISPR_Guide_Capture"]) if (config["modality_flags"]["CRISPR_Guide_Capture"]!="")and(len(config["vis_features"]["CRISPR_Guide_Capture"])>0) else None, caption="../report/dot_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
        dot_plots_Custom = report(expand(os.path.join(config["result_path"],'{{split}}','counts','plots','NORMALIZED_dot_plot_{category}_{modality}.png'),
                             category=config["vis_categories"],
                             modality=config["modality_flags"]["Custom"]) if (config["modality_flags"]["Custom"]!="")and(len(config["vis_features"]["Custom"])>0) else None, caption="../report/dot_plot.rst", category="visualization_{}".format(config["project_name"]), subcategory="{split}"),
    resources:
        mem=config.get("mem", "16G"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","dot_plot_{split}.log"),
    params:
        partition=config.get("partition"),
        ab_flag = config["modality_flags"]['Antibody_Capture'],
        crispr_flag = config["modality_flags"]['CRISPR_Guide_Capture'],
        custom_flag = config["modality_flags"]['Custom'],
        vis_categories = config["vis_categories"],
        vis_gene_lists = config["vis_gene_lists"],
        ab_features = config["vis_features"]['Antibody_Capture'],
        crispr_features = config["vis_features"]['CRISPR_Guide_Capture'],
        custom_features = config["vis_features"]['Custom'],
    script:
        "../scripts/dot_plot.R"