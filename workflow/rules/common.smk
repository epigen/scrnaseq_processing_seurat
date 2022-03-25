##### utility functions #####

def get_sample_paths(wildcards):
    return [annot.loc[wildcards.sample,'path'], annot.loc[wildcards.sample,'metadata']]
