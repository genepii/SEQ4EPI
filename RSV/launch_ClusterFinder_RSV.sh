bash ClusterFinder_RSV.sh \
    /path/to/input.fasta \  # Input FASTA file
    /path/to/output_prefix \  # Prefix for output files
    15400 \  # Fragment size
    0.00007 \  # Clustering threshold
    /path/to/reference.fasta \  # Reference FASTA file nextclade
    /path/to/genome_annotation.gff3 \  # Reference Genome annotation file (GFF3)
    /path/to/metadata.csv \  # Metadata CSV file
    /path/to/output.treefile \  # Output file for the phylogenetic tree
    /path/to/output_visualisation_folder # Folder for visualization results