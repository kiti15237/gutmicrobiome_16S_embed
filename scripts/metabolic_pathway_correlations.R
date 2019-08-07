library(vegan)
library(KEGGREST)
library(gtools)
library(RColorBrewer)
library(cowplot)
library(pheatmap)

setwd("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces")
data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"
pathway_dir = paste(data_dir, "/pathways/", sep = "")

#################################
### Read in sequence ids file ###
#################################
seqs_fasta <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/piph/seqs_filter.07_piph.fasta", quote="\"", comment.char="")
headers <- gsub(">", "", as.character(seqs_fasta[seq(1, nrow(seqs_fasta), 2), ]))
seqs <- as.character(seqs_fasta[seq(2, nrow(seqs_fasta), 2), ])
seq_id_df <- data.frame(headers)
rownames(seq_id_df) <- seqs

v <- readRDS("pca_embedded_taxa.rds")
rownames(v) <- as.character(seq_id_df[rownames(v), ])


####################################
###  Read qual vecs  ###############
####################################

glove_emb <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/AG_new/embed/glove_emb_AG_newfilter.07_100.txt",
                        quote="\"", comment.char="", row.names = 1, sep = " ", header = F)

glove_emb = glove_emb[-which(rownames(glove_emb) == '<unk>'), ]
rownames(glove_emb) <- as.character(seq_id_df[rownames(glove_emb), ])


###################################
###   Get pathway table   #########
###################################

pathway_table <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/AG_new/pathways/otu_pathway_table.RDS")
keep <- colSums(pathway_table) > 0
keep2 <- colSums(pathway_table) < nrow(pathway_table)
pathway_table <- pathway_table[, keep & keep2]


#####################################
### Match, clean, and center   ######
#####################################

embed_table_glove <- glove_emb
colnames(embed_table_glove) <- paste("dim", seq(1, ncol(embed_table_glove)), sep = "")
embed_table_pca <- v[,1:100]
colnames(embed_table_pca) <- paste("pca", seq(1, ncol(embed_table_pca)), sep = "")


taxa_names <- intersect(rownames(pathway_table), rownames(embed_table_glove)) #should be the same regardless of the embedding table
pathway_table <- pathway_table[taxa_names, ]
embed_table_glove <- embed_table_glove[taxa_names, ]
embed_table_pca <- embed_table_pca[taxa_names, ]
embed_table_pca <- apply(embed_table_pca, 2, function(x) return((x - mean(x) ) / sd(x)))
embed_table_glove <- apply(embed_table_glove, 2, function(x) return((x - mean(x) ) / sd(x)))



########################################################################
### Find correlations between embedding dimensions and pathway table####
########################################################################

getCorMat <- function(embedding_table, pathway_table){
  cor_list <- list()
  for(i in seq(1, ncol(embedding_table))){
    cor = apply(pathway_table, 2, function(pathway_vec) return(cor(pathway_vec, embedding_table[ ,i])))
    cor_list[[i]] <- cor
  }
  cor_mat <- data.frame(matrix(unlist(cor_list), byrow = T, nrow = length(cor_list)))
  colnames(cor_mat) <- colnames(pathway_table)
  rownames(cor_mat) <- colnames(embedding_table)
  return(cor_mat)
}

makeNullPathwayTable_reorder <- function(pathway_table){
  new_order <- sample(seq(1, nrow(pathway_table)), size = nrow(pathway_table))
  return(pathway_table[new_order, ])
}



#############################################################
######################## Plot heatmaps  ######################
##############################################################

null_pathway_table <- makeNullPathwayTable_reorder(pathway_table)
cor_mat_pca <- getCorMat(embed_table_pca, pathway_table)
cor_mat_glove <- getCorMat(embed_table_glove, pathway_table)
cor_mat_pca_null <- getCorMat(embed_table_pca, null_pathway_table)
cor_mat_glove_null <- getCorMat(embed_table_glove, null_pathway_table)

breaksList = seq(-0.35, 0.35, by = .01)
colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(length(breaksList))
labs<- c("pca", "embed", "pca_null", "embed_null")
i = 1
heatmaps <- list()
for(mat in list(cor_mat_pca,cor_mat_glove,cor_mat_pca_null, cor_mat_glove_null  )){
  #f = paste("C:/Users/ctata/Documents/Lab/quality_vectors/figures/", labs[i], "_pathway_heatmap.pdf", sep = "" )
  #print(f)
  #pdf(f, width = 20, height = 15)
  if(i == 2){
    legend = T
  }else{
    legend = F
  }
  if(i ==1  | i == 2 | i ==3 | i == 4){
    col_labels = rep(" ", ncol(mat))
  }
  row_labels <- rep(" ", nrow(mat))
  
  heatmaps[[i]] <- pheatmap(mat, color = colors, breaks = breaksList,
                            treeheight_col = 0, treeheight_row = 0, 
                            labels_row = row_labels, labels_col = col_labels, 
                            border_color = NA, legend = F)
  #dev.off()
  i = i + 1
}

p <- plot_grid(heatmaps[[1]][[4]], heatmaps[[2]][[4]],
               heatmaps[[3]][[4]], heatmaps[[4]][[4]])
p

pdf("../../../figures/cor_metabolic_pathways_grid.pdf", width = 5, height = 5)
p
dev.off()

setwd("C:/Users/ctata/Documents/Lab/quality_vectors/figures")
pdf("legend.pdf", width = 5, height = 5)
tmp <- pheatmap(mat, color = colors, breaks = breaksList,
         treeheight_col = 0, treeheight_row = 0, 
         labels_row = row_labels, labels_col = col_labels, 
         border_color = NA, legend = T)
tmp
dev.off()

title <- ggdraw() + draw_label("Embedding dimensions correlate with metabolic pathways", fontface='bold')
plot_grid(title, p, ncol=1, rel_heights=c(0.1, .3)) # rel_heights values control title margins
#heatmap.2(as.matrix(rbind(cor_mat_glove_df, cor_mat_pca_df)))
pheatmap(as.matrix(cor_mat_glove))


makeNullPathwayTable <- function(pathway_table){
  ps <- colMeans(pathway_table)
  nullCols <- lapply(ps, function(p) return(rbinom(n = nrow(pathway_table), size = 1, p = p)))
  nullPathwayTable <- data.frame(matrix(unlist(nullCols), byrow = F, ncol = length(nullCols)))
  colnames(nullPathwayTable) <- names(nullCols)
  return(nullPathwayTable)
}


getMaxCorr <- function(embed_vec, pathway_table){
  set.seed(123)
  #Find pathway with the highest correlation
  corrs <- apply(pathway_table, 2, cor, embed_vec)
  max_corr <- max(corrs)
  max_inx <- which(corrs == max_corr)
  
  #If we do the exact same process with randomly generated data multiple times, what are the chances we see a correlation has 
  #high as we did?
  maxPerm = 10000
  null_max_corrs <- c()
  for(iter in seq(1, maxPerm)){
    nullTable <- makeNullPathwayTable_reorder(pathway_table)
    corrs_null <- apply(nullTable, 2, cor, embed_vec)
    max_corr_null <- max(abs(corrs_null), na.rm = T)
    max_inx <- which(abs(corrs_null) == abs(max_corr_null))
    null_max_corrs <- c(null_max_corrs, max_corr_null)
  }
  pval <- sum(abs(null_max_corrs) >= abs(max_corr), na.rm=T) / maxPerm
  pval

  return(list(max_corr = max_corr,
              max_inx = max_inx,
              null_dist = null_max_corrs,
              pval = pval))
}


#####################################################################
###########  Calculate and save p values from permuation test #######
#####################################################################
increment <- 10
numGroups <- 100 / increment
#for(i in seq(1, numGroups)){
#  start <- (i-1) * increment + 1
#  end <- (i * increment)
#  corr_matches <- lapply(seq(start,end), function(i){
#    print(i)
#    return(getMaxCorr(embedding_table[ , i], pathway_table))
#  })
#  file <- paste("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/piph/corr_matches_pca_", start, "_", end, ".RDS", sep = "")
#  print(file)
#  saveRDS(corr_matches, file)
#}



####################################################################################
#################   load corr_matches, which were calculated in increments, ########
#################    merge, and resave as one object                       #########
####################################################################################

corr_matches <- list()
for(i in seq(1, numGroups)){
  start <- (i-1) * increment + 1
  end <- (i * increment)
  file <- paste("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/AG_new/pathways/corr_matches_", start, "_", end, ".RDS", sep = "")
  print(file)
  corr_matches_tmp <- readRDS(file)
  corr_matches <- c(corr_matches, corr_matches_tmp)
}

saveRDS(corr_matches, paste(pathway_dir, "corr_matches.rds", sep = ""))
saveRDS(pathway_table, paste(pathway_dir, "pathway_table.RDS", sep = ""))


##########################################################################
#########   Load objects to get stats without recalculating  #############
##########################################################################
corr_matches_glove <- readRDS(paste(pathway_dir, "corr_matches_glove.rds", sep = ""))
corr_matches_pca <- readRDS(paste(pathway_dir, "corr_matches_pca.rds", sep = ""))
unlist(lapply(corr_matches_glove, function(x) return(x$pval)))
unlist(lapply(corr_matches_pca, function(x) return(x$pval))) 

#From this exercise, I've convinced myself that there must be some real correspondence between the dimensions in embedding space and 
#some real biological function
#Get the id of the pathway each embedding dimension aligned to

pathway_id_match <- lapply(corr_matches, function(x) return(x$max_inx)) #names are embedding dimension
names(pathway_id_match) <- paste("embed_", seq(1, length(corr_matches)), sep = "")


#which pathway got picked a lot?
pathway_hits <- table(names(unlist(lapply(corr_matches, function(x) return(x$max_inx)))))
#out of 148 possible biological pathways, 66 had high correspondence with an embedding dimension or more. 
#There are 78 Desert pathways that are almost always present or almost always absent

#Get pathway id to name dictionary
lentries <- lapply(paste("map", names(pathway_hits), sep = ""), function(path) return(keggGet(path)))
pathway_names <- unlist(lapply(lentries, function(entry) return(entry[[1]]$NAME)))
pathway_id_name_dict <- as.list(pathway_names)
names(pathway_id_name_dict) <- names(pathway_hits)

#Get names of all the pathways that match an embedding dimension
embedding_pathways <- lapply(pathway_id_match, function(x) return(pathway_id_name_dict[names(x)[1]]))
names(embedding_pathways) <- names(pathway_id_match)


#Save properties and their corresponding metabolic pathways
df <- data.frame(dim = names(embedding_pathways), 
                 pathway_id =sapply(strsplit(names(unlist(embedding_pathways)), split = "\\."), `[`, 2) ,
                 pathway_names = unlist(embedding_pathways))

write.table(df, "piph/dim_pathway_dict.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#Get distance matrix
  # dot product
  cosine_sim <- cosine(t(glove_emb))
  #magnitudes
  mags = apply(glove_emb, 1, function(vec) return(sqrt(sum(vec^2))))
  
  tmp = apply(dot_product, 1, function(vec) return(vec / mags))
  cosine_sim = apply(tmp, 2, function(vec) return(vec / mags))
  
  cosine_dist = 1 - cosine_sim
  
  cosine_dist <- as.dist(cosine_dist)
  
#Permanova
perm = adonis(cosine_dist ~ Phylum, tax_table, permutation = 5)
  
embed <- function(otu, qual_vecs){
  return(otu %*% qual_vecs)
}
  

##CCA
#Samples in metadata space
#Samples in property spac

#otu table
otu = read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/otu_filtered_AG_07perc_feces.csv", 
                 sep = "\t", row.names = 1, header = T)

otu_use = otu[rownames(glove_emb), ]
otu_use = t(otu_use)
rownames(otu_use) <- gsub("X", "", rownames(otu_use))




#mapping 
mapping = read.delim2("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/AG_mapping.txt", stringsAsFactors=FALSE, row.names = 1)
mapping_use <- mapping[data.frame(mapping)$HMP_SITE == "FECAL", ]
sample_names <- intersect(rownames(mapping_use) , rownames(otu_use))
mapping_use <- mapping_use[sample_names, ]
mapping_use <- mapping_use[ , !(grepl("VIOSCREEN", colnames(mapping_use)))]
mapping_use <- data.frame(mapping_use)



otu_use <- otu_use[sample_names, ]

#Embed
embeded_otu <- embed(as.matrix(otu_use), as.matrix(glove_emb))


#check
rownames(embeded_otu) == rownames(mapping_use)

#formula
vars = c("IBD", "EXERCISE_FREQUENCY", "SEX", "ONE_LITER_OF_WATER_A_DAY_FREQUENCY", 
         "SEAFOOD_FREQUENCY", "PROBIOTIC_FREQUENCY", "OLIVE_OIL", "FRUIT_FREQUENCY", 
         "SLEEP_DURATION", "SUGAR_SWEETENED_DRINK_FREQUENCY", "MILK_CHEESE_FREQUENCY",
         "RED_MEAT_FREQUENCY","MEAT_EGGS_FREQUENCY", "VEGETABLE_FREQUENCY")

embeded_otu <- embeded_otu - min(embeded_otu)

ps <- phyloseq(otu_table(embeded_otu, taxa_are_rows = F), sample_data(mapping_use))
samples_cosine_dist <- 1 - cosine(t(otu_use))
plotCCA(ps, as.dist(samples_cosine_dist), color = "IBD", type = "species")

form <- formula(paste("embeded_otu ~ ", paste(vars, collapse = "+"), sep = ""))
cca_obj <- cca(form, data = mapping_use)

#Taxa in pathway space
#Taxa in property space









#Check out the desert
tmp <- pheatmap(cor_mat_glove)
splits <- cutree(tmp$tree_col, h =  sort(tmp$tree_col$height, decreasing = T)[7])
desert_pathways <- names(splits[splits == 2])

hist(colSums(pathway_table[ , desert_pathways]), col = "red", breaks = 50)
hist(colSums(pathway_table[ , -which(colnames(pathway_table) %in% desert_pathways)]), col = "blue", add = T, breaks = 50)






