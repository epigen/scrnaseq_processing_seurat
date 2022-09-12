# Configuration

You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. Always use absolute paths. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed
- sample annotation (sample_annotation): CSV file consisting of three columns
    -  sample_name: name of the sample (tip: keep it short, but unique)
    -  path (2 options):
        - 10X Genomics output: absolute path to the directory containing the Cell Ranger output folder filtered_feature_bc_matrix/
        - MTX files: absolute path to the directory containing the following 3 files:
            - matrix.mtx containing the counts
            - barcodes.tsv containing the cell barcodes in the first column without header (TSV)
            - features.tsv containing the feature/gene names in the first column without header (TSV)
    -  metadata (optional): abosulte path to sample metadata as CSV with the first column being cell barcodes and every other coloumn metadata for the respective barcode/cell
