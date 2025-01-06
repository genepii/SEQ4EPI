#!/bin/bash

bash ClusterFinder_FLU.sh \
    /path/to/input_HA.fasta \                      # Input FASTA file for HA sequences
    /path/to/input_NA.fasta \                      # Input FASTA file for NA sequences
    path/to/output_prefix \              		# Prefix for output files
    1700 \                                        # HA sequence length
    1400 \                                        # NA sequence length
    0.00035 \                                     # Similarity threshold for clustering HA
    0.00025 \                                     # Similarity threshold for clustering NA
    /path/to/reference_HA/ \                      # Directory for reference data (HA)
    /path/to/metadata_HA.csv \                    # Metadata CSV file for HA
    /path/to/reference_NA/ \                      # Directory for reference data (NA)
    /path/to/metadata_NA.csv \                    # Metadata CSV file for NA
    /path/to/output_HA_treefile \                  # Output file for the HA phylogenetic tree
    /path/to/output_NA_treefile \                  # Output file for the NA phylogenetic tree
    /path/to/output_visualisation_ha/ \            # Folder for HA visualisation results
    /path/to/output_visualisation_na/              # Folder for NA visualisation results