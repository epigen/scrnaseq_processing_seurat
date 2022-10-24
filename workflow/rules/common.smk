##### utility functions #####

def get_sample_paths(wildcards):
    return [annot.loc[wildcards.sample,'path'], annot.loc[wildcards.sample,'metadata']]

def get_gene_list_path(wildcards):
    return gene_list_dict[wildcards.gene_list]