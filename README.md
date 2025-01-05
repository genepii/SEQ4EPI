# SEQ4EPI

# ClusterFinder Pipeline
This repository contains a comprehensive pipeline for sequence alignment, distance calculation, clustering, and generating final tables of results. The pipeline integrates multiple tools and scripts written in Python, R, and Shell scripts.

# Table of Contents
# Prerequisites
# Installation
# Usage
# Pipeline Details
# Scripts Description
# Examples
# Contact


# Prerequisites
Ensure the following tools and packages are installed on your system:

Python 3.x: Make sure you have Python installed along with the required packages:
pandas
argparse
logging
R: Make sure you have R installed along with the required packages:
tidyverse
ggplot2
IQ-TREE: A software for phylogenetic tree construction.
Nextalign: A tool for sequence alignment.
Nextclade: A tool for sequence analysis.
TreeCluster: A tool for clustering sequences based on phylogenetic trees.

# Installation
Clone this repository to your local machine: 
       git clone https://github.com/yourusername/ClusterFinder.git
       cd ClusterFinder


# Usage
The pipeline can be executed using the provided shell script launch_treecluster.sh. Before running the script, update it with the correct paths for your input files and parameters.

To launch the pipeline: 
        ./launch_ClusterFinder.sh


Pipeline Details
Sequence Alignment with Nextalign
Sequence Analysis with Nextclade
Phylogenetic Tree Construction with IQ-TREE
Sequence Clustering with TreeCluster
Merging Metadata with Variants, Insertions, Deletions, and Clusters
Final Table Generation with R Scripts
Visualization Preparation for ITOL


# Scripts Description
Python Scripts
ClusterFinder.py: The main script for running the sequence alignment, analysis, and clustering steps.

# Usage: 

ClusterFinder.py [-h] --input_fasta INPUT_FASTA --output_prefix OUTPUT_PREFIX --genome_length GENOME_LENGTH --threshold THRESHOLD --reference REFERENCE --annotation ANNOTATION --metadata_file METADATA_FILE


--input_fasta: Path to the input FASTA file containing unaligned sequences.
--output_prefix: Prefix for the output files.
--genome_length: Length of the reference genome.
--threshold: Threshold for clustering.
--reference: Reference genome file for Nextalign.
--annotation: Annotation GFF3 file for Nextalign.
--metadata_file: Metadata CSV file.


R Scripts

final_table_of_clusters.R: Script for generating the final table of clusters.

# Usage:
Rscript final_table_of_clusters.R input_metadata.csv output_final_table.csv
  
  # input_metadata.csv: Path to the merged metadata CSV file.
  # output_final_table.csv: Path to the output final table CSV file.


Itol.R: Script for preparing the visualization data for ITOL.

# Usage:
Rscript Itol.R input_treefile.tree output_final_table.csv output_visualization_folder


 # input_treefile.tree: Path to the input tree file in Newick format.
 # output_final_table.csv: Path to the output final table CSV file.
 # output_visualization_folder: Path to the folder for storing visualization files.


 Shell Scripts
ClusterFinder.sh: Shell script for running the clustering process.

# Usage:
./ClusterFinder.sh input_fasta prefix genome_length threshold reference gff3 metadata final_table treefile output_visualisation

 # input_fasta: Path to the input FASTA file.
 # prefix: Prefix for the output files.
 # genome_length: Length of the reference genome.
 # threshold: Threshold for clustering.
 # reference: Path to the reference genome FASTA file.
 # gff3: Path to the GFF3 annotation file.
 # metadata: Path to the metadata CSV file.
 # final_table: Path to the final table CSV file.
 # treefile: Path to the output tree file in Newick format.
 # output_visualisation: Path to the folder for storing visualization files.


 launch_ClusterFinder.sh: Script to launch the entire pipeline.
# Usage:
./launch_treecluster.sh

Examples
Running ClusterFinder Pipeline
Ensure you have updated launch_treecluster.sh with the correct paths and parameters. An example of the updated script might look like this:

./Treecluster.sh \
  /path/to/input.fasta \
  /path/to/output_prefix \
  29903 \
  0.000085 \
  /path/to/reference.fasta \
  /path/to/genome_annotation.gff3 \
  /path/to/metadata.csv \
  /path/to/final_table_merged_metadata.csv \
  /path/to/output.treefile \
  /path/to/output_visualisation_folder


After updating the paths, run the script: ./launch_ClusterFinder.sh


Please ensure you have all the required dependencies installed and paths correctly set up in the scripts before running the pipeline. This README provides an overview and usage guide for the ClusterFinder pipeline. For more detailed instructions and troubleshooting, refer to the comments within each script.