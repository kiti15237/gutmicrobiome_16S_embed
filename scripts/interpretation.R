library(vegan)

#Read qual vecs
glove_emb <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/glove_emb_AG_new07perc_feces_250.txt",
                        quote="\"", comment.char="", row.names = 1, sep = " ", header = F)

glove_emb = glove_emb[-which(rownames(glove_emb) == '<unk>'), ]


#read tax table
tax_table = read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/greengenes/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt",
                       row.names = 1, header = T, sep = "\t")

tax_table_list = apply(tax_table, 1, function(tax_name) return(strsplit(as.character(tax_name), ';')))
tax_table_mat = matrix(unlist(tax_table_list), ncol = 7, byrow=T)
rownames(tax_table_mat) <- rownames(tax_table)
colnames(tax_table_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


#line them up
rownames(tax_table_mat)
tax_table  <- tax_table_mat[rownames(glove_emb), ]
tax_table = data.frame(tax_table)


#Get pathway table
pathway_table <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/allbodysites/piphillan/otu_pathway_table.RDS")

taxa_names <- intersect(rownames(pathway_table), rownames(glove_emb))
pathway_table <- pathway_table[taxa_names, ]
embedding_table <- glove_emb[taxa_names, ]

keep <- colSums(pathway_table) > 0
keep2 <- colSums(pathway_table) < nrow(pathway_table)
pathway_table <- pathway_table[, keep & keep2]


#Find correlations between embedding dimensions and pathway table
cor_list <- list()
for(i in seq(1, ncol(embedding_table))){
  cor = apply(pathway_table, 2, function(pathway_vec) return(cor(pathway_vec, embedding_table[ ,i])))
  cor_list[[i]] <- cor
}


makeNullPathwayTable <- function(pathway_table){
  ps <- colMeans(pathway_table)
  nullCols <- lapply(ps, function(p) return(rbinom(n = nrow(pathway_table), size = 1, p = p)))
  nullPathwayTable <- data.frame(matrix(unlist(nullCols), byrow = F, ncol = length(nullCols)))
  colnames(nullPathwayTable) <- names(nullCols)
  return(nullPathwayTable)
}

makeNullPathwayTable_reorder <- function(pathway_table){
  new_order <- sample(seq(1, nrow(pathway_table)), size = nrow(pathway_table))
  return(pathway_table[new_order, ])
}

getMaxCorr <- function(embed_vec, pathway_table){
  #Find pathway with the highest correlation
  corrs <- apply(pathway_table, 2, cor, embed_vec)
  max_corr <- max(corrs)
  max_inx <- which(corrs == max_corr)
  
  #If we do the exact same process with randomly generated data multiple times, what are the chances we see a correlation has 
  #high as we did?
  maxPerm = 1000
  null_max_corrs <- c()
  for(iter in seq(1, maxPerm)){
    nullTable <- makeNullPathwayTable_reorder(pathway_table)
    corrs_null <- apply(nullTable, 2, cor, embed_vec)
    max_corr_null <- max(corrs_null, na.rm = T)
    max_inx <- which(corrs_null == max_corr_null)
    null_max_corrs <- c(null_max_corrs, max_corr_null)
  }
  pval <- sum(null_max_corrs >= max_corr, na.rm=T) / maxPerm
  pval

  return(list(max_corr = max_corr,
              max_inx = max_inx,
              null_dist = null_max_corrs,
              pval = pval))
}

corr_matches <- lapply(seq(1,250), function(i){
  print(i)
  return(getMaxCorr(embedding_table[ , i], pathway_table))
})

saveRDS(corr_matches, "C:/Users/ctata/Documents/Lab/quality_vectors/interpretation/allbodysites/pathways/corr_matches.RDS")
saveRDS(pathway_table, "C:/Users/ctata/Documents/Lab/quality_vectors/interpretation/allbodysites/pathways/pathway_table.RDS")
saveRDS(embedding_table, "C:/Users/ctata/Documents/Lab/quality_vectors/interpretation/allbodysites/pathways/embedding_table.RDS")
unlist(lapply(corr_matches, function(x) return(x$pval)))

#From this exercise, I've convinced myself that there must be some real correspondence between the dimensions in embedding space and 
#some real biological function
#Get the id of the pathway each embedding dimension aligned to
pathway_id_match <- lapply(corr_matches, function(x) return(x$max_inx))


#which pathway got picked a lot?
pathway_hits <- table(names(unlist(lapply(corr_matches, function(x) return(x$max_inx)))))
#out of 137 possible biological pathways, 78 had high correspondence with an embedding dimension or more

#Get pathway id to name dictionary
lentries <- lapply(paste("map", names(pathway_hits), sep = ""), function(path) return(keggGet(path)))
pathway_names <- unlist(lapply(lentries, function(entry) return(entry[[1]]$NAME)))
pathway_id_name_dict <- as.list(pathway_names)
names(pathway_id_name_dict) <- names(pathway_hits)

#Get names of all the pathways that match an embedding dimension
embedding_pathways <- lapply(pathway_id_match, function(x) return(pathway_id_name_dict[names(x)[1]]))


df <- data.frame(topic = paste("topic_", seq(1,250), sep = ""), 
                 pathway_id = names(unlist(embedding_pathways)), pathway_names = unlist(embedding_pathways))

write.table(df, "topic_pathway_dict.txt", sep = "\t", quote = F, row.names = F, col.names = T)

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











