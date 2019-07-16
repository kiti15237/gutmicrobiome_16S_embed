
#from KEGG_genomes.R
taxa_by_ko <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/allbodysites/taxa_by_ko.rds")
otu_genome <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/allbodysites/piphillan/200p_otu_genome_hit_table.txt", header = T)
otu_genome$OTU <- gsub("seq", "", otu_genome$OTU)
otu_genome_sorted <- otu_genome[match(rownames(taxa_by_ko), otu_genome$NNgenome), ]
otu_genome_sorted$NNgenome == rownames(taxa_by_ko)
rownames(taxa_by_ko) <- otu_genome_sorted$OTU


glove_emb <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/glove_emb_AG_new07perc_feces_250.txt",
                        quote="/", comment.char="", row.names = 1, sep = " ", header = F)
glove_emb = glove_emb[-which(rownames(glove_emb) == '<unk>'), ]


taxa_names <- intersect(rownames(glove_emb), rownames(taxa_by_ko))
print(length(taxa_names))

qual_vecs <- glove_emb[taxa_names, ]
taxa_by_ko <- taxa_by_ko[taxa_names, ]
taxa_by_ko <- taxa_by_ko[ , colSums(taxa_by_ko) > 0]

matches <- c()
corr_vec <- c()
for(i in seq(1,250)){
  corrs <- cor(qual_vecs[,i], taxa_by_ko, method = "pearson")
  max_corr <- max(corrs, na.rm = T)
  inx <- which(corrs == max_corr)
  matches <- c(matches, inx)
  corr_vec <- c(corr_vec, max_corr)
}


mapGeneToPathway <- function(organism) {
  KEGG_PATHWAY_LINK_BASE <- "http://rest.kegg.jp/link/pathway/"
  pathway_link_REST_url <- paste(KEGG_PATHWAY_LINK_BASE, organism, sep="")
  
  gene_pathway <- data.frame()
  
  for (line in readLines(pathway_link_REST_url)) {
    tmp <- strsplit(line, "\t")[[1]]
    gene <- tmp[1]
    gene <- strsplit(gene, ":")[[1]][2]  
    pathway_id<- strsplit(tmp[2], organism)[[1]][2]
    
    if (is.null(gene_pathway[gene, 1])) {
      gene_pathway[gene,1] = pathway_id
    } else {
      if (is.na(gene_pathway[gene,1])) {
        gene_pathway[gene,1] = pathway_id
      } else {
        gene_pathway[gene,1] = paste(gene_pathway[gene, 1], pathway_id, sep=";")
      }
    }
  }
  names(gene_pathway) <- "pathway_id"
  gene_pathway
}


mapGeneToPathway(genomes[1])





