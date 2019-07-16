library(ape)
library(stylo)
tree <- ape::read.tree(file = "C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/97_otus.tree")
labels <- unlist(lapply(strsplit(tree$tip.label, "[.]"), function(x) return(x[1])))
tree$tip.label = labels


#otu_filtered <- read.delim("C:/Users/ctata/Documents/Lab/quality_vectors/data/silva/otu_filtered.csv", row.names = 1)
#taxa = as.character(otu_filtered$X.OTU.ID)


glove_emb <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/qual_vecs_feces.csv",
                                      quote="\"", comment.char="", row.names = 1, sep = ",")
#taxa = rownames(glove_emb_distalgut_100)
#taxa = taxa[!(taxa %in% c("CP004008", "CP004008", "CP004008" ,"CP004008","CP003344", "AC145156"))]
#taxa = taxa[taxa != "<unk>"]

#taxa = scan("taxa.txt", sep = '\t', what = character())
taxa = rownames(glove_emb)

delete = tree$tip.label[!(tree$tip.label %in% taxa)]
delete_dup = tree$tip.label[duplicated(tree$tip.label)]

tree_small <- drop.tip(tree, c(delete, delete_dup)) #How did the tree get rooted?
write.tree(tree_small, "C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/97_otus_pruned.tre")

phy_dists <- cophenetic.phylo(tree_small)
write.table(as.matrix(phy_dists), "C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/phy_dists.txt")




write(rownames(phy_dists), )
phy_dists_vec <- phy_dists[lower.tri(phy_dists, diag = F)]


glove_emb_distalgut_100_prune <- glove_emb_distalgut_100[taxa, ]
cosine_dists <- dist.cosine(as.matrix(glove_emb_distalgut_100_prune))

end = 1000000
phy_dists_small = phy_dists_vec[1:end]
cosine_dists_small = cosine_dists[1:end]
plot(phy_dists_small, cosine_dists_small )
line = lm(cosine_dists_small ~ phy_dists_small, 
          data = as.data.frame(cbind(phy_dists_small, cosine_dists_small)))
abline(line, col = "red")



cor.test(phy_dists_small, cosine_dists_small)
#We observe a moderate but robust correlation (p=4*10^-16)
#showing that evolutionary distance is coarsely related to dissimilarity between preferred neighborhoods 
#Cosine distance is a way of quantifying preference in neighboorhood. Taxa with similar vectors (small cosine distance) prefer similar neighborhoods, and are likely to be found in similar company.
#If we assume that taxonomy is a coarse measure of function (ie. same species have similar core genomes to a small extent), then we would expect to see a small but 
#robust correlation between neighborhood preference similarity and taxonomic relatedness.