import csv
import numpy as np
import random

f = open('/nfs3/PHARM/David_Lab/christine/public_data/AG_new/dada2_output/filtered_150/seqtab_final_filter.07_99agg.txt', 'r')
outfile = open('/nfs3/PHARM/David_Lab/christine/public_data/AG_new/glove/glove_input_filter.07_99agg.txt', mode = 'w')
test_samples_file = open('/nfs3/PHARM/David_Lab/christine/public_data/AG_new/dada2_output/filtered_150/test_samples_filter.07.txt', 'w')
print("filter .07")


writer = csv.writer(outfile, delimiter = "\t", quoting = csv.QUOTE_NONE, escapechar = '')

taxa_names = f.readline()
taxa_names = taxa_names.strip().split("\t")
i = 0
test_samples = []
random.seed(0)
for line in f:
    vec = line.split("\t")
    sample_id = vec[0]
    if random.random() > 0.15:	
        present = [float(i) > 0 for i in vec[1:]]
        writer.writerow(np.array(taxa_names)[present])
        print(i, end = '\t')
        i = i + 1
    else:
        test_samples.append(sample_id)
	
test_samples_file.write("\t".join(test_samples))
f.close()
outfile.close()
