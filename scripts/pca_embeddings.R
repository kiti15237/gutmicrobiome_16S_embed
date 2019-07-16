library(text2vec)
library(ggfortify)

qual_vecs <- read.table("C:/Users/ctata/Documents/Lab/Microbiome_ASD_16S/Microbiome_ASD_16S/embeddings/glove_emb_AG_newfilter.07_100.txt",
                        quote="\"", comment.char="", row.names = 1)
qual_vecs <- qual_vecs[rownames(qual_vecs) != "<unk>", ]

pca <- prcomp(qual_vecs)


kpcof <- read.csv("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/annotations/kpcof_annotation.txt",
                        sep="", colClasses = rep("character", 6))
kpcof[is.na(kpcof)] = ""


annotations <- read.csv("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/annotations/genus_species_annotation.txt",
                        sep="", colClasses = rep("character", 3))
annotations[is.na(annotations)] = ""

pairs <- data.frame(Genus1 = c("Bacteroides", "Roseburia"), Species1 = c("uniformis", "intestinalis"),
                    Genus2 = c("Escherichia/Shigella", "Blautia"), Species2 = c("coli", "hydrogenotrophica"),
                    stringsAsFactors = F)
                    
seqs_pairs1 <- lapply(pairs$Species1, function(x) return(annotations[annotations$Species == x , "Seq"]))

seqs_pairs2 <- lapply(pairs$Species2, function(x) return(annotations[annotations$Species == x , "Seq"]))

crossfeeding_dists <- mapply(function(seqs1, seqs2){
  mat1 <- as.matrix(qual_vecs[seqs1, ])
  mat2 <-  as.matrix(qual_vecs[seqs2, ])
  dists <- 1 - sim2( mat1, mat2, method = "cosine")
  return(dists)
}, seqs_pairs1, seqs_pairs2)


rowWise_cosine_dist <- function(mat){
  mat <- as.matrix(mat)
  dot <- mat %*% t(mat)
  mags <- apply(mat, 1, function(x) return(sqrt(sum(x^2))))
  dot1 <- apply(dot, 1, function(x) return(x/mags))
  dot2 <- apply(dot1, 2, function(x) return(x/mags))
  return(1 - dot2)
}

out_dists <- rowWise_cosine_dist(qual_vecs)

qual_vecs$size <- rep(1, nrow(qual_vecs))
qual_vecs$alpha <- rep(0.3, nrow(qual_vecs))
qual_vecs$shape <- rep(21, nrow(qual_vecs))
for(i in seq(1, nrow(pairs))){
  qual_vecs$label[rownames(qual_vecs) %in% seqs_pairs1[[i]]] <- paste(pairs$Genus1[i], pairs$Species1[i])
  qual_vecs$label[rownames(qual_vecs) %in% seqs_pairs2[[i]]] <- paste(pairs$Genus2[i], pairs$Species2[i])
  qual_vecs$size[rownames(qual_vecs) %in% seqs_pairs1[[i]]] <- 2
  qual_vecs$size[rownames(qual_vecs) %in% seqs_pairs2[[i]]] <- 2
  
  qual_vecs$alpha[rownames(qual_vecs) %in% seqs_pairs1[[i]]] <- 1
  qual_vecs$alpha[rownames(qual_vecs) %in% seqs_pairs2[[i]]] <- 1
  
  qual_vecs$shape[rownames(qual_vecs) %in% seqs_pairs1[[i]]] <- "pair1"
  qual_vecs$shape[rownames(qual_vecs) %in% seqs_pairs2[[i]]] <- "pair2"
}

autoplot(pca, data = qual_vecs, colour = 'label', size = 'size', alpha = 'alpha', shape = 'shape')


qual_vecs <- read.table("C:/Users/ctata/Documents/Lab/Microbiome_ASD_16S/Microbiome_ASD_16S/embeddings/glove_emb_AG_newfilter.07_250.txt",
                        quote="\"", comment.char="", row.names = 1)
qual_vecs <- qual_vecs[rownames(qual_vecs) != "<unk>", ]
qual_vecs <- qual_vecs[order(rownames(qual_vecs)), ]
kpcof <- kpcof[match(rownames(qual_vecs), kpcof$Seq), ]
annotations <- annotations[match(rownames(qual_vecs), annotations$Seq), ]
sum(kpcof$Seq == rownames(qual_vecs))

pca <- prcomp(qual_vecs)
tmp <- qual_vecs
tmp$Phylum <- kpcof$Phylum
tmp$Class <- kpcof$Class
tmp$Genus <- annotations$Genus
tmp$color <- rep("NA", nrow(tmp))

mark <- c("Bacteroides", "Roseburia", "Blautia", "Clostridium", "")
markGenus <- annotations$Genus[annotations$Genus %in% mark]
tmp$color[annotations$Genus %in% mark] <- markGenus
tmp$alpha <- rep(0.1, nrow(tmp))
tmp$alpha[annotations$Genus %in% mark] <- 1
tmp$size <- rep(1, nrow(tmp))
tmp$size[annotations$Genus %in% mark] <- 3
tmp$shape <- rep(1, nrow(tmp))
tmp$shape[annotations$Species == ""]


autoplot(pca, data = tmp, colour = 'color', alpha = 'alpha', size = 'size')
