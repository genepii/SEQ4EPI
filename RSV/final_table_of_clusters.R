argv <- commandArgs(TRUE)

input_table <- as.character(argv[1])
output_table <- as.character(argv[2])

table <- read.csv(input_table, colClasses = "character")

clusttable <- aggregate(table$cluster, by = list(table$cluster, table$location, table$deletions), function(x) unique(x))

colnames(clusttable) <- c("cluster", "location", "deletions", "FinalCluster")

subcluster <- c(LETTERS, letters)

for (currcluster in unique(clusttable$cluster)) {
  clusttable[clusttable$cluster == currcluster, ]$FinalCluster <- paste0(clusttable[clusttable$cluster == currcluster, ]$FinalCluster, subcluster[1:nrow(clusttable[clusttable$cluster == currcluster, ])])
}

table <- merge(table, clusttable, by = c("cluster", "location", "deletions"))

table$FinalCluster <- ifelse(table$cluster == "-1", "NA", table$FinalCluster)

write.table(table, output_table, row.names = F, sep = ";", quote = F)
