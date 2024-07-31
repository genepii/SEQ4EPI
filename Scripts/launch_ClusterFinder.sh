

./ClusterFinder.sh \
  /home/.../file.fasta \ #file fasta 
  /home/.../output_Folder \ #results folder
  29903 \ #length of reference genome
 0.000085 \ #threshold for the assignment clusters
  /home/.../nextclade_dataset/reference.fasta \ #reference genome fasta
  /home/.../nextclade_dataset/genome_annotation.gff3 \ #reference genome gff for query sequences alignment with nextclade
  /home/.../metadata.csv \ #metadata 
  /home/.../final_table_merged_metadata.csv \ #final table of results with clusters assigned
  /home/.../output.treefile \ #output tree in newick format
  /home/.../output_visualisation_Folder #visualisation folder