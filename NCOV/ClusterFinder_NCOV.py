import os
import subprocess
import pandas as pd
import argparse
import logging
import glob
from Bio import SeqIO

# Set up logging
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
            return
        else:
            logging.info(f"MAFFT completed successfully: {result.stdout.decode()}")
    except Exception as e:
        logging.error(f"Error during MAFFT execution: {str(e)}")
        return

# Function to filter sequences that are too short, have too much ambiguous data, or too many mixed sites
def filter_sequences(input_fasta, output_fasta, log_file, min_length=26000, max_missing_data=4000, max_mixed_sites=20):
    logging.info(f"Filtering sequences in {input_fasta}, generating {output_fasta}")
    filtered_sequences = []

    # Bases representing mixed sites
    mixed_sites = set(['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'])

    sequence_count = 0
    valid_sequence_count = 0
    short_sequence_count = 0
    too_many_missing_data_count = 0
    too_many_mixed_sites_count = 0

    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence_count += 1
            seq = str(record.seq)

            # Count missing data and mixed sites
            missing_data_count = seq.count("N") + seq.count("-")
            mixed_site_count = sum([1 for base in seq if base in mixed_sites])

            logging.debug(f"Analyzing sequence {record.id}: length = {len(seq)}, "
                          f"missing data = {missing_data_count}, mixed sites = {mixed_site_count}")

            # Apply filtering criteria
            if len(seq) < min_length:
                short_sequence_count += 1
                logging.debug(f"Sequence {record.id} rejected: too short (length = {len(seq)}).")
            elif missing_data_count > max_missing_data:
                too_many_missing_data_count += 1
                logging.debug(f"Sequence {record.id} rejected: too much missing data ({missing_data_count}).")
            elif mixed_site_count > max_mixed_sites:
                too_many_mixed_sites_count += 1
                logging.debug(f"Sequence {record.id} rejected: too many mixed sites ({mixed_site_count}).")
            else:
                # Write valid sequence to the filtered FASTA file
                SeqIO.write(record, output_handle, "fasta")
                valid_sequence_count += 1

            # Add to filtered sequences list for potential report
            filtered_sequences.append({
                'seqName': record.id,
                'sequence_length': len(seq),
                'missing_data_count': missing_data_count,
                'mixed_site_count': mixed_site_count
            })

    # Logs after filtering
    logging.info(f"Total number of sequences analyzed: {sequence_count}")
    logging.info(f"Number of valid sequences written to {output_fasta}: {valid_sequence_count}")
    logging.info(f"Number of sequences too short (< {min_length} bases): {short_sequence_count}")
    logging.info(f"Number of sequences with too much missing data (> {max_missing_data}): {too_many_missing_data_count}")
    logging.info(f"Number of sequences with too many mixed sites (> {max_mixed_sites}): {too_many_mixed_sites_count}")

    # Log filtered sequences to a CSV file for analysis
    log_filtered_sequences(filtered_sequences, log_file)
    logging.info(f"Filtered FASTA file generated: {output_fasta}")

# Function to log filtered sequences to a CSV file
def log_filtered_sequences(filtered_sequences, log_file):
    if filtered_sequences:
        logging.info(f"Recording filtered sequences to {log_file}...")
        df = pd.DataFrame(filtered_sequences)
        df.to_csv(log_file, index=False)
    else:
        logging.info("No sequences filtered.")

def run_nextalign(input_fasta, output_dir, reference, annotation):
    check_tool("nextalign")
    logging.info("Aligning sequences with Nextalign")
    nextalign_command = f"nextalign run -r {reference} -g {annotation} -O {output_dir} {input_fasta}"
    subprocess.run(nextalign_command, shell=True, check=True)
    logging.info(f"Nextalign finished. Results in: {output_dir}")

    # List generated files
    logging.info(f"Files generated in {output_dir}: {os.listdir(output_dir)}")

def run_nextclade(input_fasta, output_dir):
    check_tool("nextclade")
    logging.info("Analyzing sequences with Nextclade...")
    nextclade_dataset = "/path/to/nextclade_dataset"  # Update with actual path to Nextclade dataset
    nextclade_command = f"nextclade run -D {nextclade_dataset} -O {output_dir} {input_fasta}"
    result = subprocess.run(nextclade_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logging.error(f"Nextclade failed with message: {result.stderr.decode()}")
    else:
        logging.info(f"Nextclade finished. Results in: {output_dir}")
        logging.info(f"Standard output of Nextclade: {result.stdout.decode()}")
        logging.info(f"Nextclade errors: {result.stderr.decode()}")
        logging.info(f"Files generated in {output_dir}: {os.listdir(output_dir)}")

def build_phylogenetic_tree(aligned_fasta, output_prefix):
    logging.info("Building phylogenetic tree with IQ-TREE...")
    iqtree_checkpoint = f"{output_prefix}.ckp.gz"
    
    if os.path.exists(iqtree_checkpoint):
        os.remove(iqtree_checkpoint)

    iqtree_command = f"iqtree -s {aligned_fasta} -m GTR+G -bb 1000 -nt AUTO --seed 12345 --keep-ident -pre {output_prefix} -redo"

    try:
        result = subprocess.run(iqtree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logging.error(f"IQ-TREE failed: {result.stderr.decode()}")
            return
        else:
            logging.info(f"IQ-TREE completed successfully: {result.stdout.decode()}")
    except Exception as e:
        logging.error(f"Error during IQ-TREE execution: {str(e)}")
        return

    treefile_pattern = f"{output_prefix}.treefile"
    treefile_list = glob.glob(treefile_pattern)

    if not treefile_list:
        logging.error(f"No .treefile found with prefix {output_prefix}.")
        return
    
    treefile = treefile_list[0]

    if os.path.exists(treefile):
        logging.info(f"Phylogenetic tree file found: {treefile}")
    else:
        logging.error(f"The tree file {treefile} was not generated.")
        return

def run_treecluster(tree_file, aligned_fasta, output_file, threshold):
    check_tool("python")
    treecluster_script = "/path/to/TreeCluster.py"
    logging.info("Clustering sequences with TreeCluster...")
    treecluster_command = f"python {treecluster_script} -i {tree_file} -o {output_file} -t {threshold}"
    try:
        subprocess.run(treecluster_command, shell=True, check=True)
        logging.info(f"Clustering finished. Results file: {output_file}")
    except Exception as e:
        logging.error(f"Error during TreeCluster execution: {str(e)}")

def merge_metadata_with_variants(metadata_file, nextclade_csv_file, nextalign_csv_file, cluster_file, output_file):
    logging.info("Merging metadata with Nextclade, Nextalign, and cluster results...")

    # Check if files exist before continuing
    if not os.path.exists(metadata_file):
        logging.error(f"The file {metadata_file} does not exist.")
        return
    if not os.path.exists(nextclade_csv_file):
        logging.error(f"The file {nextclade_csv_file} does not exist.")
        return
    if not os.path.exists(nextalign_csv_file):
        logging.error(f"The Nextalign file {nextalign_csv_file} does not exist.")
        return
    if not os.path.exists(cluster_file):
        logging.error(f"The file {cluster_file} does not exist.")
        return

    metadata = pd.read_csv(metadata_file, encoding='utf-8-sig')
    logging.info(f"Columns in metadata: {metadata.columns.tolist()}")

    try:
        nextclade_data = pd.read_csv(nextclade_csv_file, delimiter=';', on_bad_lines='skip')
        logging.info(f"Columns in nextclade_data: {nextclade_data.columns.tolist()}")
    except pd.errors.ParserError as e:
        logging.error(f"Error reading Nextclade CSV file: {e}")
        return

    try:
        nextalign_data = pd.read_csv(nextalign_csv_file, delimiter=',', on_bad_lines='skip')
        logging.info(f"Columns in nextalign_data: {nextalign_data.columns.tolist()}")
    except pd.errors.ParserError as e:
        logging.error(f"Error reading Nextalign CSV file: {e}")
        return

    clusters = pd.read_csv(cluster_file, sep='\t', names=["SequenceName", "cluster"])
    clusters.rename(columns={"SequenceName": "seqName"}, inplace=True)
    logging.info(f"Columns in clusters: {clusters.columns.tolist()}")
    
    metadata['seqName'] = metadata['seqName'].astype(str)
    nextclade_data['seqName'] = nextclade_data['seqName'].astype(str)
    nextalign_data['seqName'] = nextalign_data['seqName'].astype(str)
    clusters['seqName'] = clusters['seqName'].astype(str)

    # Merge data
    combined_data = pd.merge(metadata, nextclade_data, on='seqName', how='outer', suffixes=('', '_clade'))
    combined_data = pd.merge(combined_data, nextalign_data, on='seqName', how='outer', suffixes=('', '_align'))
    combined_data = pd.merge(combined_data, clusters, on='seqName', how='outer')

    required_columns = ['seqName', 'collection_date', 'location', 'deletions', 'insertions', 'cluster', 'clade_display', 'clade_who', 'clade_nextstrain', 'partiallyAliased', 'Nextclade_pango', 'qc.overallScore']
    missing_columns = [col for col in required_columns if col not in combined_data.columns]
    if missing_columns:
        logging.error(f"Missing columns after merge: {missing_columns}")
        return

    combined_data = combined_data[required_columns]

    combined_data['deletions'] = combined_data['deletions'].fillna('')
    combined_data['insertions'] = combined_data['insertions'].fillna('')
    combined_data['cluster'] = combined_data['cluster'].fillna(-1)

    combined_data['cluster_group'] = combined_data.groupby(['cluster', 'location', 'deletions', 'insertions']).ngroup()

    def assign_alphabetical_labels(group):
        group = group.copy()
        group['Final Cluster'] = [f"{group['cluster'].iloc[0]}{chr(65 + i)}" for i in range(len(group))]
        return group

    combined_data = combined_data.groupby(['cluster', 'location', 'deletions', 'insertions'], group_keys=False).apply(assign_alphabetical_labels)

    combined_data.to_csv(output_file, index=False)
    logging.info(f"Merged files saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Complete pipeline for alignment, filtering, clustering, and generating the final table.')
    parser.add_argument('--input_fasta', required=True, help='FASTA file containing the unaligned sequences.')
    parser.add_argument('--output_prefix', required=True, help='Prefix for output files.')
    parser.add_argument('--genome_length', type=int, required=True, help='Genome length.')
    parser.add_argument('--threshold', type=float, required=True, help='Clustering threshold.')
    parser.add_argument('--metadata_file', required=True, help='CSV file containing metadata.')
    parser.add_argument('--reference', required=True, help='Path to the reference file (e.g., reference.fasta).')
    parser.add_argument('--input_annotation', required=True, help='Path to the annotation file (e.g., genome_annotation.gff3).')

    args = parser.parse_args()

    fasta_basename = os.path.splitext(os.path.basename(args.input_fasta))[0]

    mafft_output_aligned_fasta = f"{args.output_prefix}_{fasta_basename}_mafft.aligned.fasta"
    filtered_fasta = f"{args.output_prefix}_{fasta_basename}_filtered.aligned.fasta"
    filtered_log_file = f"{args.output_prefix}_{fasta_basename}_filtered_sequences_log.csv"
    nextalign_output_dir = f"{args.output_prefix}_{fasta_basename}_nextalign"
    nextclade_output_dir = f"{args.output_prefix}_{fasta_basename}_nextclade"
    nextalign_aligned_fasta = os.path.join(nextalign_output_dir, 'nextalign.aligned.fasta')
    nextalign_insertions_csv = os.path.join(nextalign_output_dir, 'nextalign.insertions.csv')
    merged_metadata_file = f"{args.output_prefix}_merged_metadata.csv"
    iqtree_output_prefix = f"{args.output_prefix}_{fasta_basename}_iqtree"

    run_mafft(args.input_fasta, mafft_output_aligned_fasta)
    filter_sequences(mafft_output_aligned_fasta, filtered_fasta, filtered_log_file)
    run_nextalign(filtered_fasta, nextalign_output_dir, args.reference, args.input_annotation)
    run_nextclade(filtered_fasta, nextclade_output_dir)
    build_phylogenetic_tree(mafft_output_aligned_fasta, iqtree_output_prefix)
    run_treecluster(iqtree_output_prefix + ".treefile", mafft_output_aligned_fasta, f"{args.output_prefix}_clusters.tsv", args.threshold)
    merge_metadata_with_variants(args.metadata_file, nextclade_insertions_csv, nextalign_insertions_csv, f"{args.output_prefix}_clusters.tsv", merged_metadata_file)

if __name__ == "__main__":
    main()
