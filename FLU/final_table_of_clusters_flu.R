# Récupération des arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Lecture du fichier d'entrée
# Ajustez 'sep' selon votre fichier d'entrée ("," ou ";")
table <- read.csv(input_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Colonnes définissant le pattern interne, incluant cluster_NA
pattern_cols <- c("location_x", "cluster_HA", "deletions_x", "insertions_x", 
                  "cluster_NA", "deletions_y", "insertions_y")

# Colonnes finales attendues
expected_columns <- c("seqName", "collection_date_x", "location_x", "index_x",
                      "clade_HA", "subclade", "qc.overallScore_x", "qc.overallStatus_x", 
                      "cluster_HA", "totalSubstitutions_HA", "clade_NA", "totalSubstitutions_NA", 
                      "deletions_x", "insertions_x", "cluster_NA", "deletions_y", "insertions_y")

# Vérification des colonnes nécessaires
all_needed_cols <- unique(c(expected_columns, pattern_cols))
missing_columns <- setdiff(all_needed_cols, colnames(table))
if (length(missing_columns) > 0) {
  stop("Colonnes manquantes dans la table : ", paste(missing_columns, collapse = ", "))
}

# Créer la colonne FinalCluster
table$FinalCluster <- NA

# Lettres pour distinguer les différents patterns
subcluster_letters <- c(LETTERS, letters)

# Ajouter la colonne 'Different' pour capturer les patterns empêchant l'attribution
table$Different <- NA

# Parcourir chaque cluster_HA
for (curr_cluster in unique(table$cluster_HA)) {
  # Lignes correspondant à ce cluster_HA
  cluster_rows <- which(table$cluster_HA == curr_cluster)
  
  # Extraire les données du pattern
  cluster_data <- table[cluster_rows, pattern_cols, drop = FALSE]
  
  # Créer une clé unique par ligne pour identifier les patterns
  patterns <- apply(cluster_data, 1, function(x) paste(x, collapse = "|"))
  unique_patterns <- unique(patterns)
  
  if (length(unique_patterns) == 1) {
    # Un seul pattern dans ce cluster_HA
    table$FinalCluster[cluster_rows] <- paste0(curr_cluster, "A")
  } else {
    # Plusieurs patterns, attribuer une lettre par pattern
    for (p_idx in seq_along(unique_patterns)) {
      pattern_rows <- cluster_rows[patterns == unique_patterns[p_idx]]
      letter <- subcluster_letters[p_idx]
      table$FinalCluster[pattern_rows] <- paste0(curr_cluster, letter)
    }
  }
}

# Identifier les séquences avec NA dans FinalCluster et déterminer les patterns différents
table$Different[is.na(table$FinalCluster)] <- apply(
  table[is.na(table$FinalCluster), pattern_cols, drop = FALSE], 1, function(row) {
    differences <- sapply(pattern_cols, function(col) {
      if (!is.na(row[col])) {
        paste0(col, ": ", row[col])
      } else {
        NA
      }
    })
    paste(na.omit(differences), collapse = " | ")
  }
)
# Remplacer les clusters avec moins de 3 occurrences par NA
cluster_counts <- table(table$FinalCluster)
clusters_to_na <- names(cluster_counts[cluster_counts < 3])
table$FinalCluster[table$FinalCluster %in% clusters_to_na] <- NA

# Remplacer toutes les occurrences de -1 et NA dans FinalCluster par NA
table$FinalCluster[table$FinalCluster == "-1" | table$FinalCluster == "NA"] <- NA

# S'assurer que toutes les colonnes attendues sont présentes
for (col in expected_columns) {
  if (!col %in% colnames(table)) {
    table[[col]] <- NA
  }
}

# Conserver uniquement les colonnes demandées + FinalCluster
final_cols <- c(expected_columns, "FinalCluster")
final_table <- table[, final_cols, drop = FALSE]

# Écriture de la table finale
write.table(final_table, output_file, row.names = FALSE, sep = ";", quote = FALSE)
