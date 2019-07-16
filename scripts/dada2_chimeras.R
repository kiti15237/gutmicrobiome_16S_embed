library(dada2); packageVersion("dada2")

seqtabs = list()
for(i in seq(40)){
	if(i < 10){
		i = paste("0", i, sep = "")
	}
	file = paste("../dada2_output/filtered_150/seqtab_0", i, ".rds",  sep = "")
	print(file)
	seqtabs[[i]] <- readRDS(file)

}

all <- mergeSequenceTables(tables = seqtabs)

seqtab <- removeBimeraDenovo(all, method = "consensus", multithread = TRUE)

saveRDS(seqtab, "../dada2_output/filtered_150/seqtab_final.rds")
