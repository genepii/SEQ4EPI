import os
import subprocess
import pandas as pd
import argparse
import logging
from Bio import SeqIO

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def check_tool(tool_name):
    result = subprocess.run(f"which {tool_name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logging.error(f"{tool_name} is not installed or is not in the PATH.")
        raise FileNotFoundError(f"{tool_name} is not installed or is not in the PATH.")
    else:
        logging.info(f"{tool_name} found: {result.stdout.decode().strip()}")

def run_mafft(input_fasta, output_aligned_fasta, tmpdir="/path/to/tmp"):
    check_tool("mafft")
    logging.info("Aligning sequences with MAFFT...")
    mafft_command = f"mafft --auto {input_fasta} > {output_aligned_fasta}"
    env = os.environ.copy()
    env["TMPDIR"] = tmpdir
    try:
        result = subprocess.run(mafft_command, shell=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logging.error(f"MAFFT failed: {result.stderr.decode()}")
        else:
            logging.info("MAFFT completed successfully.")
    except Exception as e:
        logging.error(f"Error executing MAFFT: {str(e)}")

def run_nextclade(input_fasta, output_dir, dataset_path):
    check_tool("nextclade")
    logging.info(f"Running Nextclade analysis for dataset {dataset_path}...")
    nextclade_command = f"nextclade run -D {dataset_path} -O {output_dir} {input_fasta}"
    result = subprocess.run(nextclade_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logging.error(f"Nextclade failed: {result.stderr.decode()}")
    else:
        logging.info(f"Nextclade finished for {dataset_path}. Results are in: {output_dir}")

def build_phylogenetic_tree(aligned_fasta, output_prefix):
    logging.info("Building phylogenetic tree with IQ-TREE...")
    iqtree_command = f"iqtree -s {aligned_fasta} -m GTR+G -bb 1000 -nt AUTO --seed 12345 --keep-ident -pre {output_prefix} -redo"
    try:
        result = subprocess.run(iqtree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logging.error(f"IQ-TREE failed: {result.stderr.decode()}")
        else:
            logging.info(f"IQ-TREE completed successfully.")
    except Exception as e:
        logging.error(f"Error executing IQ-TREE: {str(e)}")

def run_treecluster(tree_file, output_file, threshold):
    check_tool("python")
    treecluster_script = "/path/to/TreeCluster.py"
    logging.info("Clustering sequences with TreeCluster...")
    treecluster_command = f"python {treecluster_script} -i {tree_file} -o {output_file} -t {threshold}"
    try:
        subprocess.run(treecluster_command, shell=True, check=True)
        logging.info(f"Clustering completed. Result file: {output_file}")
    except Exception as e:
        logging.error(f"Error executing TreeCluster: {str(e)}")

def filter_sequences(input_fasta, output_fasta, log_file, min_length=1400, max_missing_data=50, max_mixed_sites=10):
    logging.info(f"Filtering sequences in {input_fasta}...")
    filtered_sequences = []
    mixed_sites = set(['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'])
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            missing_data_count = seq.count("N") + seq.count("-")
            mixed_site_count = sum(1 for base in seq if base in mixed_sites)
            if len(seq) >= min_length and missing_data_count <= max_missing_data and mixed_site_count <= max_mixed_sites:
                SeqIO.write(record, output_handle, "fasta")
                filtered_sequences.append({
                    'seqName': record.id,
                    'sequence_length': len(seq),
                    'missing_data_count': missing_data_count,
                    'mixed_site_count': mixed_site_count
                })
    pd.DataFrame(filtered_sequences).to_csv(log_file, index=False)
    logging.info(f"Filtered FASTA file generated: {output_fasta}")

def merge_metadata_with_variants(metadata_file_ha, metadata_file_na, nextclade_csv_file_ha, nextclade_csv_file_na, cluster_file_ha, cluster_file_na, output_file):
    logging.info("Merging metadata with Nextclade HA, Nextclade NA results, and clusters...")
    metadata_ha = pd.read_csv(metadata_file_ha, dtype={'seqName': str})
    metadata_na = pd.read_csv(metadata_file_na, dtype={'seqName': str})
    nextclade_data_ha = pd.read_csv(nextclade_csv_file_ha, delimiter=';', dtype={'seqName': str}).rename(columns={
        'clade': 'clade_HA', 'totalSubstitutions': 'totalSubstitutions_HA'})
    nextclade_data_na = pd.read_csv(nextclade_csv_file_na, delimiter=';', dtype={'seqName': str}).rename(columns={
        'clade': 'clade_NA', 'totalSubstitutions': 'totalSubstitutions_NA'})
    clusters_ha = pd.read_csv(cluster_file_ha, sep='\t', names=["seqName", "cluster_HA"], dtype={'seqName': str})
    clusters_na = pd.read_csv(cluster_file_na, sep='\t', names=["seqName", "cluster_NA"], dtype={'seqName': str})

    combined_data = metadata_ha.merge(nextclade_data_ha, on='seqName', how='outer') \
                               .merge(clusters_ha, on='seqName', how='outer') \
                               .merge(metadata_na, on='seqName', how='outer') \
                               .merge(nextclade_data_na, on='seqName', how='outer') \
                               .merge(clusters_na, on='seqName', how='outer')
    combined_data.to_csv(output_file, index=False)
    logging.info(f"Merged files saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Analysis pipeline for HA and NA.")
    
    # Arguments for HA and NA fragment FASTA files
    parser.add_argument("--input_fasta_ha", required=True, help="FASTA file for HA fragment.")
    parser.add_argument("--input_fasta_na", required=True, help="FASTA file for NA fragment.")
    
    # Arguments for the genome length of each fragment
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files.")
    parser.add_argument('--genome_length_ha', type=int, required=True, help='Genome length for HA.')
    parser.add_argument('--genome_length_na', type=int, required=True, help='Genome length for NA.')
    
    # Arguments for clustering thresholds for each fragment
    parser.add_argument('--threshold_ha', type=float, required=True, help='Threshold for HA clustering.')
    parser.add_argument('--threshold_na', type=float, required=True, help='Threshold for NA clustering.')
    
    # Arguments for metadata files for each fragment
    parser.add_argument('--metadata_file_ha', required=True, help='CSV file containing metadata for HA.')
    parser.add_argument('--metadata_file_na', required=True, help='CSV file containing metadata for NA.')
    
    # Arguments for Nextclade datasets for each fragment
    parser.add_argument('--dataset_ha', required=True, help='Nextclade dataset for HA.')
    parser.add_argument('--dataset_na', required=True, help='Nextclade dataset for NA.')

    # Parse the arguments
    args = parser.parse_args()

    # Define output directories and file names
    output_dir_ha = f"{args.output_prefix}_HA"
    output_dir_na = f"{args.output_prefix}_NA"
    filtered_fasta_ha = f"{output_dir_ha}_filtered.fasta"
    filtered_fasta_na = f"{output_dir_na}_filtered.fasta"
    aligned_fasta_ha = f"{output_dir_ha}_aligned.fasta"
    aligned_fasta_na = f"{output_dir_na}_aligned.fasta"
    tree_file_ha = f"{output_dir_ha}.treefile"
    tree_file_na = f"{output_dir_na}.treefile"
    cluster_file_ha = f"{output_dir_ha}_clusters.txt"
    cluster_file_na = f"{output_dir_na}_clusters.txt"
    nextclade_output_ha = os.path.join(output_dir_ha, 'nextclade.csv')
    nextclade_output_na = os.path.join(output_dir_na, 'nextclade.csv')
    merged_metadata_file = f"{args.output_prefix}_merged_metadata.csv"

    # Processing steps for HA fragment
    filter_sequences(args.input_fasta_ha, filtered_fasta_ha, f"{output_dir_ha}_log.csv")
    run_nextclade(filtered_fasta_ha, output_dir_ha, args.dataset_ha)
    run_mafft(filtered_fasta_ha, aligned_fasta_ha)
    build_phylogenetic_tree(aligned_fasta_ha, output_dir_ha)
    run_treecluster(tree_file_ha, cluster_file_ha, args.threshold_ha)

    # Processing steps for NA fragment
    filter_sequences(args.input_fasta_na, filtered_fasta_na, f"{output_dir_na}_log.csv")
    run_nextclade(filtered_fasta_na, output_dir_na, args.dataset_na)
    run_mafft(filtered_fasta_na, aligned_fasta_na)
    build_phylogenetic_tree(aligned_fasta_na, output_dir_na)
    run_treecluster(tree_file_na, cluster_file_na, args.threshold_na)

    # Merge metadata and clustering results
    merge_metadata_with_variants(
        args.metadata_file_ha, 
        args.metadata_file_na,
        nextclade_output_ha, 
        nextclade_output_na,
        cluster_file_ha, 
        cluster_file_na,
        merged_metadata_file
    )

if __name__ == "__main__":
    main()
