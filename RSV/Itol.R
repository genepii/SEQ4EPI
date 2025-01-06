# Charger les bibliothèques nécessaires
library(ape)
library(dplyr)
library(ggtree)
library(ggplot2)
library(randomcoloR)

# Récupérer les arguments de la ligne de commande
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  stop("Usage: Rscript Visualisation.R <tree_file.treefile> <data_file.csv> <output_directory>")
}

tree_file <- args[1]
data_file <- args[2]
output_dir <- args[3]

# Vérification et création du répertoire si nécessaire
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Lire l'arbre phylogénétique au format Newick
tree <- read.tree(tree_file)

# Lire le fichier CSV contenant les informations sur les séquences et les clusters
data <- read.csv(data_file, header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE)

# Convertir seqName en caractère pour éviter les incompatibilités de type
data$seqName <- as.character(data$seqName)

# Vérifier la structure des données pour s'assurer que les colonnes nécessaires existent
if(!("seqName" %in% names(data)) | !("FinalCluster" %in% names(data))) {
  stop("Le fichier CSV doit contenir les colonnes 'seqName' et 'FinalCluster'")
}

# Extraire le nom de fichier sans extension pour le titre
tree_file_name <- tools::file_path_sans_ext(basename(tree_file))
plot_title <- paste("Maximum-likelihood phylogenies and clustering assignments of", tree_file_name)

# Associer les données de cluster aux pointes de l'arbre
tree$tip.label <- as.character(tree$tip.label)
tip_data <- data.frame(label = tree$tip.label)
tip_data <- left_join(tip_data, data, by = c("label" = "seqName"))

# Filtrer les clusters avec au moins 3 occurrences
cluster_counts <- data %>% group_by(FinalCluster) %>% tally()
clusters_with_at_least_three <- cluster_counts %>% filter(n >= 3) %>% pull(FinalCluster)

# Générer des couleurs uniquement pour les clusters ayant au moins 3 occurrences
generate_colors <- function(clusters) {
  unique_clusters <- unique(clusters)
  cluster_colors <- distinctColorPalette(length(unique_clusters))
  names(cluster_colors) <- unique_clusters
  return(cluster_colors)
}

# Générer des couleurs pour les clusters ayant au moins 3 occurrences
colors <- generate_colors(clusters_with_at_least_three)

# Attribuer NA (aucune couleur) aux clusters avec moins de 3 occurrences
tip_data$color <- sapply(tip_data$FinalCluster, function(cluster) {
  if (cluster %in% clusters_with_at_least_three) {
    return(colors[cluster])
  } else {
    return(NA)  # Pas de couleur pour les clusters avec moins de 3 occurrences
  }
})

# Créer un vecteur de formes pour les clusters (formes pleines)
shapes <- 15:(15 + length(clusters_with_at_least_three) - 1)
names(shapes) <- clusters_with_at_least_three

# Associer les formes uniquement aux clusters ayant au moins 3 occurrences
tip_data$shape <- sapply(tip_data$FinalCluster, function(cluster) {
  if (cluster %in% clusters_with_at_least_three) {
    return(shapes[cluster])
  } else {
    return(NA)  # Pas de forme pour les clusters avec moins de 3 occurrences
  }
})

# Générer un fichier PNG de l'arbre coloré avec labels et formes
output_png_file <- file.path(output_dir, "colored_tree_with_labels_and_shapes.png")
p <- ggtree(tree) %<+% tip_data +
  geom_tiplab(aes(label = label), show.legend = FALSE, size = 3) +  # Ajuster la taille du texte si nécessaire
  geom_tippoint(aes(color = FinalCluster, shape = FinalCluster), size = 3, na.rm = TRUE) +  # na.rm = TRUE pour ne pas afficher les clusters sans couleur ni forme
  scale_color_manual(values = colors, na.translate = FALSE) +  # na.translate = FALSE pour ne pas inclure NA dans la légende
  scale_shape_manual(values = shapes, na.translate = FALSE) +
  theme_tree2() +
  ggtitle(plot_title) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title = "Clusters"), shape = guide_legend(title = "Clusters"))

# Ajuster la taille de l'image
ggsave(output_png_file, plot = p, width = 20, height = 15, units = "in")  # Ajuster les dimensions si nécessaire

# Exporter l'arbre avec les couleurs pour iTOL
output_tree_file <- file.path(output_dir, "colored_tree_for_iTOL.nwk")
write.tree(tree, file = output_tree_file)

# Créer un fichier annotation pour iTOL uniquement pour les clusters avec au moins 3 occurrences
annotation_file <- file.path(output_dir, "itol_annotations.txt")
annotation_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR SPACE",
  "DATASET_LABEL Cluster Colors",
  "COLOR #ff0000",
  "LEGEND_TITLE Cluster",
  paste("LEGEND_SHAPES", paste(rep(1, length(clusters_with_at_least_three)), collapse = " ")),
  paste("LEGEND_COLORS", paste(colors[clusters_with_at_least_three], collapse = " ")),
  paste("LEGEND_LABELS", paste(clusters_with_at_least_three, collapse = " ")),
  "STRIP_WIDTH 25",
  "MARGIN 0",
  "BORDER_WIDTH 1",
  "BORDER_COLOR #0000FF",
  "SHOW_INTERNAL 0",
  "DATA"
)

for(i in 1:nrow(data)) {
  if (data$FinalCluster[i] %in% clusters_with_at_least_three) {
    annotation_lines <- c(annotation_lines, paste(data$seqName[i], colors[data$FinalCluster[i]], "1", sep = " "))
  }
}

writeLines(annotation_lines, con = annotation_file)

cat("Arbre coloré avec labels, formes, titre et légende, fichier d'annotation, et image PNG générés avec succès dans le répertoire :", output_dir, "\n")
