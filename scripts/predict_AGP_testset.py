

import pandas as pd
import numpy as np
import helper_functions as hf
import pickle
import matplotlib.pyplot as plt
from sklearn import preprocessing
import importlib
importlib.reload(hf)
import math
import copy


data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"
fig_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/figures/"


def getQualVecs():
    qual_vec_file = data_dir + "embed/embed_.07_100dim.txt"
    qual_vecs = pd.read_csv(qual_vec_file, sep = " ", index_col = 0, dtype = {0:str})
    qual_repseqs = pd.read_csv(data_dir + "embed/seqs_.07.fasta", sep = "\t", header = None)
    
    import re
    ids = qual_repseqs.iloc[range(0, qual_repseqs.shape[0], 2), 0]
    ids = [re.sub(">", "", i) for i in ids.values]

    seqs = qual_repseqs.iloc[range(1, qual_repseqs.shape[0], 2), 0]
    seqs = [str(i) for i in seqs.values]

    #Drop <unk> character
    ids = ids[0: len(ids)-1]
    seqs = seqs[0: len(seqs)-1]
    qual_vecs = qual_vecs.iloc[0: len(seqs), :]

    print(len(ids))
    print(len(seqs))
    print(qual_vecs.shape)
    return(qual_vecs, ids, seqs)



#Load data
filt = ".07"
f = open(data_dir + "/AG_new/otu_train_" + str(filt) + ".obj", "rb")
otu_train = pickle.load(f)
f.close()

f = open(data_dir + "/AG_new/otu_test_" + str(filt) + ".obj", "rb")
otu_test = pickle.load(f)
f.close()

f = open(data_dir + "/AG_new/map_train_" + str(filt) + ".obj", "rb")
map_train = pickle.load(f)
f.close()

f = open(data_dir + "/AG_new/map_test_" + str(filt) +  ".obj", "rb")
map_test = pickle.load(f)
f.close()

qual_vecs, embed_ids, embed_seqs = getQualVecs()

otu_train = hf.matchOtuQual(otu_train, embed_ids, embed_seqs)
otu_test = hf.matchOtuQual(otu_test, embed_ids, embed_seqs)


#Normalize with asinh
target = "IBD"
f = plt.figure(figsize=(15,5))
X_train, X_val, X_test, y_train, y_val, y_test = hf.getMlInput(otu_train, otu_test, map_train, map_test, 
                                                            target = target, asinNormalized = True)
X_train = pd.concat([X_train, X_val], axis = 0)
y_train = y_train + y_val
plt.subplot(1, 2, 1)
auc_asin, auc_train_asin, fpr_asin, tpr_asin, prec_asin, f1_asin, _ = hf.predictIBD(X_train, y_train, X_test, y_test, graphTitle = "Normalized asinh Taxa Abundances " + str(X_train.shape[1]) + " features",
              max_depth = 10, n_estimators = 65, weight = 5, plot = True, plot_pr = True)

f.savefig(fig_dir + "curves_AGP_test_asin.pdf")


# In[ ]:


#Embed
f = plt.figure(figsize=(15,5))
X_train, X_val, X_test, y_train, y_val, y_test = hf.getMlInput(otu_train, otu_test, map_train, map_test, 
                                                            target = target, embed = True, qual_vecs = qual_vecs)
X_train = pd.concat([X_train, X_val], axis = 0)
y_train = y_train + y_val
plt.subplot(1, 2, 1)
auc_embed, auc_train_embed, fpr_embed, tpr_embed, prec_embed, f1_embed, _ = hf.predictIBD(X_train, y_train, X_test, y_test, graphTitle = "Embedding weighted by averaging taxa "+ str(X_train.shape[1]) + " features",
              max_depth = 10, n_estimators = 170,  weight = 5, plot = True, plot_pr = True)

f.savefig(fig_dir + "curves_AGP_test_embed.pdf")


# In[ ]:


#PCA
f = plt.figure(figsize=(15,5))
X_train, X_val, X_test, y_train, y_val, y_test = hf.getMlInput(otu_train, otu_test, map_train, map_test, 
                                                            target = target, pca_reduced = True, numComponents = 100)
X_train = pd.concat([X_train, X_val], axis = 0)
y_train = y_train + y_val
plt.subplot(1, 2, 1)
auc_pca, auc_train_pca, fpr_pca, tpr_pca, prec_pca, f1_pca, _  = hf.predictIBD(X_train, y_train, X_test, y_test, graphTitle = "PCA dimensionality reduced " + str(X_train.shape[1]) + " features", 
              max_depth = 3, n_estimators = 110, weight = 5, plot = True, plot_pr = True)
f.savefig(fig_dir + "curves_AGP_test_pca.pdf")

