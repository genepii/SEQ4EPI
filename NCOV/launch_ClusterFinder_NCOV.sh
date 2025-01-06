bash ClusterFinder_NCOV.sh \
  /path/to/input.fasta \  # Input FASTA file
  /path/to/output_prefix \  # Prefix for output files
  29903 \  # Fragment size
  0.000085 \  # Clustering size
  /path/to/reference.fasta \  # Reference FASTA file nextclade
  /path/to/genome_annotation.gff3 \  #Reference Genome annotation file (GFF3)
  /path/to/metadata.csv \  # Metadata CSV file
  /path/to/final_table_merged_metadata.csv \  # Merging Metadata with Variants, Insertions, Deletions, and Clusters
  /path/to/output.treefile \  # Output file for the phylogenetic tree
  /path/to/output_visualisation_folder # Folder for visualization results