#!/bin/bash 

# Define variables passed as arguments
fasta_ha=$1
fasta_na=$2
prefix=$3
genlength_ha=$4
genlength_na=$5
threshold_ha=$6
threshold_na=$7
reference_ha=$8
metadata_ha=$9
reference_na=${10}
metadata_na=${11}
treefile_ha=${12}
treefile_na=${13}
output_visualisation_ha=${14}
output_visualisation_na=${15}

# Ensure the viewing directory exists
echo "Variables:"
echo "Fragment HA fasta : $fasta_ha"
echo "Fragment NA fasta : $fasta_na"
echo "Prefix : $prefix"
echo "Longueur génome HA : $genlength_ha"
echo "Longueur génome NA : $genlength_na"
echo "Seuil de clustering HA : $threshold_ha"
echo "Seuil de clustering NA : $threshold_na"
echo "Référence HA : $reference_ha"
echo "Metadata HA : $metadata_ha"
echo "Référence NA : $reference_na"
echo "Metadata NA : $metadata_na"
echo "Treefile HA : $treefile_ha"
echo "Treefile NA : $treefile_na"
echo "Visualisation HA : $output_visualisation_ha"
echo "Visualisation NA : $output_visualisation_na"

# Set output folder with prefix checking
output_dir="./${prefix}_results"
mkdir -p "$output_dir"

# Check and create the view folders for HA and NA
mkdir -p "$output_visualisation_ha" || { echo "Erreur : Impossible de créer le dossier $output_visualisation_ha"; exit 1; }
mkdir -p "$output_visualisation_na" || { echo "Erreur : Impossible de créer le dossier $output_visualisation_na"; exit 1; }

python ClusterFinder_FLU.py \
    --input_fasta_ha "$fasta_ha" \
    --input_fasta_na "$fasta_na" \
    --output_prefix "${output_dir}/${prefix}" \
    --genome_length_ha "$genlength_ha" \
    --genome_length_na "$genlength_na" \
    --threshold_ha "$threshold_ha" \
    --threshold_na "$threshold_na" \
    --metadata_file_ha "$metadata_ha" \
    --metadata_file_na "$metadata_na" \
    --dataset_ha "$reference_ha" \
    --dataset_na "$reference_na"

# Generate the final cluster table for each shard
merged_metadata_file="${output_dir}/${prefix}_merged_metadata.csv"
final_table="${output_dir}/${prefix}_final_table.csv"

# Running the R script to generate the final cluster table
Rscript /path/to/final_table_of_clusters_flu.R "$merged_metadata_file" "$final_table"

# Verifying the creation of the final_table.csv file
if [[ ! -f "$final_table" ]]; then
  echo "Erreur : le tableau final des clusters n'a pas été créé."
  exit 1
fi

# Visualizing Clusters with iTOL for HA Tree
if [[ -f "$treefile_ha" ]]; then
  Rscript /path/to/Itol.R "$treefile_ha" "$final_table" "$output_visualisation_ha"
else
  echo "Erreur : le fichier arbre HA ($treefile_ha) est manquant."
fi

# Visualizing Clusters with iTOL for NA Tree
if [[ -f "$treefile_na" ]]; then
  Rscript /path/to/Itol.R "$treefile_na" "$final_table" "$output_visualisation_na"
else
  echo "Erreur : le fichier arbre NA ($treefile_na) est manquant."
fi
