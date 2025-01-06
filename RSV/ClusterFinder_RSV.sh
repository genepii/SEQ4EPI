#!/bin/bash

# Define variables passed as arguments
fasta=$1
prefix=$2
genlength=$3
threshold=$4
reference=$5
gff3=$6
metadata=$7
final_table=${prefix}_final_table.csv
treefile=$8
output_visualisation=${9}

echo "variables:"
echo $fasta
echo $prefix
echo $genlength
echo $threshold
echo $reference
echo $gff3
echo $metadata
echo $final_table
echo $treefile
echo $output_visualisation

# Ensure the viewing directory exists
if [ ! -d "$output_visualisation" ]; then
  mkdir -p "$output_visualisation"
fi

# script for clustering
python /path/to/ClusterFinder_RSV.py \
    --input_fasta "$fasta" \
    --output_prefix "$prefix" \
    --genome_length "$genlength" \
    --threshold "$threshold" \
    --reference "$reference" \
    --input_annotation "$gff3" \
    --metadata_file "$metadata"

# generate the final table
merged_metadata_file="${prefix}_merged_metadata.csv"
Rscript /path/to/final_table_of_clusters.R "$merged_metadata_file" "$final_table"

# visualisation with ITOL
Rscript /path/to/Itol.R "$treefile" "$final_table" "$output_visualisation"
