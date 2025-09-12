##### utility functions #####

def get_sample_path(wildcards):
    return annot.loc[wildcards.sample,'path']

def get_gene_list_path(wildcards):
    return gene_list_dict[wildcards.gene_list]

def get_cell_cycle_s_phase_genes(wildcards):
    if config["cell_cycle"]["s_phase_genes"]=="" or config["cell_cycle"]["s_phase_genes"]=="tirosh2015":
        return []
    else:
        return config["cell_cycle"]["s_phase_genes"]

def get_cell_cycle_g2m_phase_genes(wildcards):
    if config["cell_cycle"]["g2m_phase_genes"]=="" or config["cell_cycle"]["g2m_phase_genes"]=="tirosh2015":
        return []
    else:
        return config["cell_cycle"]["g2m_phase_genes"]


def get_module_gene_lists(wildcards):
    return [v for v in config.get("module_gene_lists", {}).values() if v != ""]


def get_vis_gene_lists(wildcards):
    return [v for v in config.get("vis_gene_lists", {}).values() if v != ""]