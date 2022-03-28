# Configuration

You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. Always use absolute paths. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed
- sample annotation (sample_annotation): CSV file consisting of three columns
    -  sample_name: name of the sample (tip: keep it short)
    -  path: absolute path to the directory containing the Cell Ranger output folder filtered_feature_bc_matrix/
    -  metadata: abosulte path to sample metadata as CSV with the first column being "barcode" (ie cells) and every other coloumn metadata for the respective barcode/cell
