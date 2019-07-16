library(KEGGREST)
library(stylo)

combineKOlists <- function(l){
  combo <- do.call(rbind, lapply(lapply(l, unlist), "[",
                        unique(unlist(c(sapply(l,names))))))
  return(combo)
}


result = tryCatch({
  expr
}, warning = function(w) {
  warning-handler-code
}, error = function(e) {
  error-handler-code
}, finally = {
  cleanup-code
}


#Read in otu to ko table
otu_genome_hit_table <- read.delim("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/allbodysites/piphillan/200p_otu_genome_hit_table.txt")
otu_genome_hit_table <- otu_genome_hit_table[!duplicated(otu_genome_hit_table$OTU), ]

genomes <- unique(otu_genome_hit_table$NNgenome)
getKOLists <- function(genomes){
  ko_lists <- list()
  
  
  #ko_lists <- ko_lists[1:365]
  #genomes <- genomes[!(genomes %in% names(ko_lists))]
  
  
  for(i in seq(along = genomes)){
    cat(paste(i, " ") )
    genome <- genomes[i]
    print(genome)
    
    errorFlag <- FALSE
    genome_entry <- tryCatch({
      keggGet(paste("genome:", genome, sep = ""))
    }, error = function(e){
      cat("error finding ")
      print(as.character(genome))
      errorFlag <<- TRUE
    })
    if(errorFlag){
      next
    }
    
    taxa <- genome_entry[[1]]$TAXONOMY
    taxa_name <- taxa$LINEAGE
    taxa_id <- taxa$TAXONOMY
    taxa_subspecies <- genome_entry[[1]]$DEFINITION
    taxa_keywords <- genome_entry[[1]]$KEYWORDS
    ko_list <- table(names(keggLink(genome_entry[[1]]$ENTRY, "ko")))
    ko_lists[[as.character(genome)]] <- ko_list
  }
  return(ko_lists)
}

#we have to do pathway presence/absence, because honestly nothing else makes sense
getPathwayTable <- function(otu_genome_hit_table){
  pathway_lists <- list()
  org_names <- c()
  for(i in seq(along= genomes)){
    print(paste(i, genomes[i]))
    
    errorFlag <- FALSE
    pathways <- tryCatch({
      unique(names(keggLink(genomes[i], "path")))
      
    }, error = function(e){
      cat("error finding ")
      print(as.character(genomes[i]))
      errorFlag <<- TRUE
    })
    if(errorFlag){
      next
    }
    pathways <- gsub(paste("path:", as.character(genomes[i]), sep = ""), "", pathways)
    df = data.frame(matrix(rep(1, length(pathways)), nrow = 1))
    
    colnames(df) <- pathways
    pathway_lists[[i]] <- df
    org_names <- c(org_names, as.character(genomes[i]))
  }
  
  pathway_table <- rbind.fill(pathway_lists)
  pathway_table[is.na(pathway_table)] <- 0
  rownames(pathway_table) <- org_names
  return(pathway_table)
}

pathway_table <- getPathwayTable(otu_genome_hit_table)
otu_genome_pruned <- otu_genome_hit_table[otu_genome_hit_table$NNgenome %in% rownames(pathway_table), ]
otu_pathway_table <- pathway_table[as.character(otu_genome_pruned$NNgenome), ]
rownames(otu_pathway_table) <- gsub("seq", "", otu_genome_pruned$OTU)

saveRDS(otu_pathway_table, "otu_pathway_table.RDS")
write.table(otu_pathway_table, "otu_pathway_table.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)






#Get ko pathway ids

taxa_by_ko <- combineKOlists(ko_lists)
taxa_by_ko[is.na(taxa_by_ko)] <- 0 #This is now a taxa by ko table, perfect for comparing to quality vectors!
saveRDS(taxa_by_ko, "C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/taxa_by_ko.rds")

taxa_by_ko <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/taxa_by_ko.rds")
#Change from genome name to otu id
otu_ids <- as.character(otu_genome_hit_table$OTU[match(rownames(taxa_by_ko), as.character(otu_genome_hit_table$NNgenome))])
#otu_ids <- otu_ids[!duplicated(otu_ids)]
#otu_ids <- otu_ids[match(rownames(taxa_by_ko), otu_ids)]


rownames(taxa_by_ko) <- otu_ids
taxa_by_ko <- taxa_by_ko[!duplicated(rownames(taxa_by_ko)), ]
rownames(taxa_by_ko) <- gsub("seq", "", rownames(taxa_by_ko))
#saveRDS(taxa_by_ko, "C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/taxa_by_ko.rds")



#Read in quality vector file
#glove_emb <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/qual_vecs.csv",
#                        quote="\"", comment.char="", row.names = 1)

qual_vecs <- read.csv("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/qual_vecs.csv", row.names=1)


#Match names between quality vectors and taxa_by_ko
qual_vecs <- qual_vecs[rownames(qual_vecs) %in% rownames(taxa_by_ko), ]
taxa_by_ko <- taxa_by_ko[rownames(taxa_by_ko) %in% rownames(qual_vecs), ]

taxa_by_ko_sorted <- taxa_by_ko[order(rownames(taxa_by_ko)), ]
qual_vecs_sorted <- qual_vecs[order(rownames(qual_vecs)), ]

rownames(taxa_by_ko_sorted) == rownames(qual_vecs_sorted)

saveRDS(taxa_by_ko_sorted, "C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/taxa_by_ko_sorted.rds")
saveRDS(qual_vecs_sorted, "C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/qual_vecs_sorted.rds")

#################################
library(stylo)
#Read in objects
taxa_by_ko_sorted <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/taxa_by_ko_sorted.rds")
qual_vecs_sorted <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/qual_vecs_sorted.rds")

#Calculate quality vector distances (cosine distances)
qual_vec_dists <- dist.cosine(as.matrix(qual_vecs_sorted))

#Calculate ko distances (cosine distances)
ko_dists <- dist.cosine(as.matrix(taxa_by_ko_sorted))


cor.test(qual_vec_dists, ko_dists)
