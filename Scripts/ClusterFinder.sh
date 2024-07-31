
fasta=$1
prefix=$2
genlength=$3
threshold=$4
reference=$5
gff3=$6
metadata=$7
final_table=$8
treefile=$9
output_visualisation=$10

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

#
#Clustering
python ClusterFinder.py --input_fasta $fasta --output_prefix $prefix --genome_length $genlength --threshold $threshold --reference $reference --annotation $gff3 --metadata_file $metadata
#R table
Rscript final_table_of_clusters.R ${prefix}_merged_metadata.csv $final_table
#Visualization_ITOL
Rscript Itol.R $treefile $final_table $output_visualisation