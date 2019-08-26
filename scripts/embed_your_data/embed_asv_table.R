library(readr)
data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/" ##CHANGE TO YOUR DATA DIRECTORY

#########################
### Read ASV table ######
#########################
asv_file = paste(data_dir, "/halfvarson/seqtab_t.txt", sep = "") ## CHANGE TO YOUR DATASET OF INTEREST
seqtab <- read.table(asv_file, row.names = 1, header = T, sep = "\t")
colnames(seqtab) <- gsub("X", "", colnames(seqtab))
seqtab <- t(seqtab)

#############################
## Read best hits table #####
## from blast ###############
#############################

best_hits <- read.delim(paste(data_dir , "/halfvarson/embed/best_hits.tsv", sep = ""), header=FALSE, row.names = 1)  ## CHANGE TO BLAST OUTPUT - SEE README
colnames(best_hits) <- c("hit_id", "query_seq", "hit_seq", "evalue", "bitscore")

#Filter best_hits table to only include hits that pass the e-value threshold
best_hits <- best_hits[best_hits$evalue < 1*10^(-29), ]

best_hits <- best_hits[as.character(best_hits$query_seq) %in% colnames(seqtab), ]

#################################
### Assign nearest neighbor id ##
#################################

#Drop any ASVs from the table that don't have near enough hits in the transformation matrix
#seqtab <- seqtab[ , colnames(seqtab) %in% rownames(best_hits)] #17784 taxa left

seqtab_hits <- seqtab[ , as.character(best_hits$query_seq) ]

#Assign the id of each ASV's nearest hit in the embedding transformation table.
colnames(seqtab_hits) <- best_hits$hit_id


##############################
## Read quality vector  ######
## transformation table ######
##############################
transform_mat_file <- paste(data_dir, "/embed/embed_.07_100dim.txt", sep = "") #CHANGE TO EMBED IN A DIFFERENT DIMENSIONAL SPACE
qual_vecs <- read.table(transform_mat_file)
qual_vecs <- qual_vecs[colnames(seqtab_hits), ]


#############################################
### Take dot product and save file ##########
#############################################

embedded <- as.matrix(asinh(seqtab_hits)) %*% as.matrix(qual_vecs)
embedded_file = paste(data_dir, "halfvarson/embed/seqtab_embedded_.07_100dim", sep = "")
saveRDS(embedded, paste(embedded_file, ".rds", sep = ""))
write.table(embedded, paste(embedded_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)








