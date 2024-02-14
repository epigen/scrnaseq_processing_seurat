##### utility functions #####

def get_sample_paths(wildcards):
    return annot.loc[wildcards.sample,'path']

def get_gene_list_path(wildcards):
    return gene_list_dict[wildcards.gene_list]

def get_large_mem_mb(wildcards, input):
    return max(2.5 * input.size_mb, int(config.get("mem", "16000")))