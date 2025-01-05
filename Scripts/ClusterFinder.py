import os
import subprocess
import pandas as pd
import argparse
import logging

# Configurer la journalisation
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_tool(tool_name):
    """ Vérifie si un outil est disponible dans le PATH """
    result = subprocess.run(f"which {tool_name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logging.error(f"{tool_name} n'est pas installé ou n'est pas dans le PATH.")
        raise FileNotFoundError(f"{tool_name} n'est pas installé ou n'est pas dans le PATH.")
    else:
        logging.info(f"{tool_name} trouvé : {result.stdout.decode().strip()}")

def run_nextalign(input_fasta, output_dir, reference, annotation):
    check_tool("nextalign")
    logging.info("Alignement des séquences avec Nextalign...")
    nextalign_command = f"nextalign run -r {reference} -g {annotation} -O {output_dir} {input_fasta}"
    subprocess.run(nextalign_command, shell=True, check=True)
    logging.info(f"Nextalign terminé. Résultats dans : {output_dir}")

    # Lister les fichiers générés
    logging.info(f"Fichiers générés dans {output_dir} : {os.listdir(output_dir)}")

def run_nextclade(input_fasta, output_dir):
    check_tool("nextclade")
    logging.info("Analyse des séquences avec Nextclade...")
    nextclade_dataset = "/home/.../nextclade_dataset"  # Mettre à jour avec le chemin réel du dataset Nextclade
    nextclade_command = f"nextclade run -D {nextclade_dataset} -O {output_dir} {input_fasta}"
    result = subprocess.run(nextclade_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logging.error(f"Nextclade a échoué avec le message : {result.stderr.decode()}")
    else:
        logging.info(f"Nextclade terminé. Résultats dans : {output_dir}")
        logging.info(f"Sortie standard de Nextclade : {result.stdout.decode()}")
        logging.info(f"Erreurs de Nextclade : {result.stderr.decode()}")
        logging.info(f"Fichiers générés dans {output_dir} : {os.listdir(output_dir)}")

def build_phylogenetic_tree(aligned_fasta, output_prefix):
    logging.info("Construction de l'arbre phylogénétique avec IQ-TREE...")
    iqtree_checkpoint = f"{output_prefix}.ckp.gz"
    if os.path.exists(iqtree_checkpoint):
        os.remove(iqtree_checkpoint)
   
    iqtree_command = f"iqtree -s {aligned_fasta} -m GTR+G -bb 1000 -nt AUTO --keep-ident -pre {output_prefix}"
    subprocess.run(iqtree_command, shell=True, check=True)
    logging.info(f"Fichier de l'arbre phylogénétique généré : {output_prefix}.treefile")

def inspect_nextclade_file(nextclade_file):
    logging.info("Inspection du fichier généré par Nextclade...")
    with open(nextclade_file, 'r') as file:
        lines = file.readlines()
        for line in lines[:5]:  # Afficher les 5 premières lignes pour inspection
            logging.info(line.strip())

def merge_metadata_with_variants(metadata_file, nextclade_csv_file, insertions_file, deletions_file, cluster_file, output_file):
    logging.info("Fusion des métadonnées avec les variants, insertions, deletions et clusters...")
    metadata = pd.read_csv(metadata_file, encoding='utf-8-sig')
    logging.info(f"Colonnes dans metadata: {metadata.columns.tolist()}")

    if not os.path.exists(nextclade_csv_file):
        logging.error(f"Le fichier {nextclade_csv_file} n'existe pas.")
        return

    inspect_nextclade_file(nextclade_csv_file)  # Inspecter le fichier généré par Nextclade

    try:
        nextclade_data = pd.read_csv(nextclade_csv_file, delimiter=';', on_bad_lines='skip')
        logging.info(f"Colonnes dans nextclade_data: {nextclade_data.columns.tolist()}")
    except pd.errors.ParserError as e:
        logging.error(f"Erreur lors de la lecture du fichier CSV : {e}")
        return

    # Convertir les colonnes 'seqName' en chaînes de caractères
    metadata['seqName'] = metadata['seqName'].astype(str)
    nextclade_data['seqName'] = nextclade_data['seqName'].astype(str)

    # Charger les insertions et les délétions
    insertions = pd.read_csv(insertions_file)
    deletions = pd.read_csv(deletions_file)
    logging.info(f"Colonnes dans insertions: {insertions.columns.tolist()}")
    logging.info(f"Colonnes dans deletions: {deletions.columns.tolist()}")
    insertions['seqName'] = insertions['seqName'].astype(str)
    deletions['seqName'] = deletions['seqName'].astype(str)

    combined_data = pd.merge(metadata, nextclade_data, left_on='seqName', right_on='seqName', how='outer', suffixes=('', '_y'))
    combined_data = pd.merge(combined_data, insertions, on='seqName', how='left', suffixes=('', '_insertions'))
    if 'deletions' in deletions.columns:
        combined_data = pd.merge(combined_data, deletions[['seqName', 'deletions']], on='seqName', how='left')
    else:
        combined_data = pd.merge(combined_data, deletions[['seqName', 'errors']], on='seqName', how='left')
        combined_data.rename(columns={'errors': 'deletions'}, inplace=True)

    # Charger les clusters générés par TreeCluster
    clusters = pd.read_csv(cluster_file, sep='\t', names=["SequenceName", "cluster"])
    clusters.rename(columns={"SequenceName": "seqName"}, inplace=True)
    logging.info(f"Colonnes dans clusters: {clusters.columns.tolist()}")
    clusters['seqName'] = clusters['seqName'].astype(str)

    combined_data = pd.merge(combined_data, clusters, on='seqName', how='outer')  # Fusionner en utilisant 'outer' pour inclure toutes les séquences

    # Vérifier les colonnes disponibles après toutes les fusions
    logging.info(f"Colonnes dans combined_data après fusion: {combined_data.columns.tolist()}")

    # Assurez-vous que toutes les colonnes nécessaires sont présentes
    required_columns = ['seqName', 'collection_date', 'location', 'deletions', 'insertions', 'cluster']
    missing_columns = [col for col in required_columns if col not in combined_data.columns]
    if missing_columns:
        logging.error(f"Colonnes manquantes après fusion: {missing_columns}")
        return

    # Sélectionner uniquement les colonnes nécessaires
    try:
        combined_data = combined_data[required_columns]
    except KeyError as e:
        logging.error(f"Erreur lors de la sélection des colonnes : {e}")
        return

    # Remplacer les valeurs non finies par une valeur par défaut
    combined_data['deletions'].fillna('', inplace=True)
    combined_data['insertions'].fillna('', inplace=True)
    combined_data['cluster'].fillna(-1, inplace=True)  # Assurez-vous que 'cluster' ne contient pas de valeurs NA

    # Ajouter une colonne pour le sous-cluster final
    combined_data['cluster_group'] = combined_data.groupby(['cluster', 'location', 'deletions', 'insertions']).ngroup()

    # Créer le Final Cluster avec les lettres consécutives
    def assign_alphabetical_labels(group):
        group = group.copy()
        group['Final Cluster'] = [f"{group['cluster'].iloc[0]}{chr(65 + i)}" for i in range(len(group))]
        return group

    combined_data = combined_data.groupby(['cluster', 'location', 'deletions', 'insertions']).apply(assign_alphabetical_labels)

    combined_data.to_csv(output_file, index=False)
    logging.info(f"Fichiers fusionnés sauvegardés : {output_file}")

def run_treecluster(tree_file, aligned_fasta, output_file, threshold):
    check_tool("python")  # Assure que python est dans le PATH pour exécuter le script TreeCluster.py
    treecluster_script = '/home/.../TreeCluster.py'  # Mett jour ce chemin avec le chemin correct
    logging.info("Clustering des séquences avec TreeCluster...")
    treecluster_command = f"python {treecluster_script} -i {tree_file} -o {output_file} -t {threshold}"
    subprocess.run(treecluster_command, shell=True, check=True)
    logging.info(f"Clustering terminé. Fichier de résultats : {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Pipeline complet pour alignement, le calcul de distances, le clustering et la génération du tableau final.')
    parser.add_argument('--input_fasta', required=True, help='Fichier FASTA contenant les séquences non alignées.')
    parser.add_argument('--output_prefix', required=True, help='Préfixe pour les fichiers de sortie.')
    parser.add_argument('--genome_length', type=int, required=True, help='Longueur du génome.')
    parser.add_argument('--threshold', type=float, required=True, help='Seuil pour le clustering.')
    parser.add_argument('--reference', required=True, help='Fichier de référence pour Nextalign.')
    parser.add_argument('--annotation', required=True, help='Fichier annotation GFF3 pour Nextalign.')
    parser.add_argument('--metadata_file', required=True, help='Fichier CSV contenant les métadonnées.')
   
    args = parser.parse_args()

    # Fichiers de sortie intermédiaires
    nextalign_output_dir = f"{args.output_prefix}_nextalign"
    aligned_fasta = os.path.join(nextalign_output_dir, "nextalign.aligned.fasta")
    insertions_file = os.path.join(nextalign_output_dir, "nextalign.insertions.csv")
    deletions_file = os.path.join(nextalign_output_dir, "nextalign.errors.csv")
    nextclade_output_dir = f"{args.output_prefix}_nextclade"
    nextclade_csv_file = os.path.join(nextclade_output_dir, "nextclade.csv")
    merged_metadata_file = f"{args.output_prefix}_merged_metadata.csv"
    tree_file = f"{args.output_prefix}.treefile"
    cluster_output_file = f"{args.output_prefix}_clusters.txt"
    final_table_file = f"{args.output_prefix}_final_table.csv"

    print(args.reference)
    print(args.annotation)
    # Étape 1 : Alignement avec Nextalign
    run_nextalign(args.input_fasta, nextalign_output_dir, args.reference, args.annotation)

    # Étape 2 : Analyse avec Nextclade
    run_nextclade(aligned_fasta, nextclade_output_dir)

    # Vérifier si le fichier CSV de Nextclade a été généré
    if not os.path.exists(nextclade_csv_file):
        logging.error(f"Le fichier {nextclade_csv_file} n'existe pas. Assurez-vous que Nextclade génère ce fichier.")
        return

    # Étape 3 : Construction de l'arbre phylogénétique avec IQ-TREE
    build_phylogenetic_tree(aligned_fasta, args.output_prefix)

    # Vérifier si le fichier de l'arbre a été généré
    if not os.path.exists(tree_file):
        logging.error(f"Le fichier {tree_file} n'a pas été généré.")
        return

    # Étape 4 : Clustering des séquences avec TreeCluster
    run_treecluster(tree_file, aligned_fasta, cluster_output_file, args.threshold)

    # Vérifier si le fichier de clustering a été généré
    if not os.path.exists(cluster_output_file):
        logging.error(f"Le fichier {cluster_output_file} n'a pas été généré.")
        return

    # Étape 5 : Fusion des métadonnées avec les variants, insertions, deletions et clusters
    merge_metadata_with_variants(args.metadata_file, nextclade_csv_file, insertions_file, deletions_file, cluster_output_file, merged_metadata_file)

    # Vérifier si le fichier fusionné a été généré
    if not os.path.exists(merged_metadata_file):
        logging.error(f"Le fichier {merged_metadata_file} n'a pas été généré.")
        return

if __name__ == "__main__":
    main()
