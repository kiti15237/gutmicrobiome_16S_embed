
# coding: utf-8

# In[1]:


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
from sklearn.model_selection import KFold

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

data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"
fig_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/figures/"

qual_vecs, embed_ids, embed_seqs = getQualVecs()



otu = pd.read_csv(data_dir + "halfvarson/seqtab.txt", sep = "\t")
otu.shape

mapping = pd.read_csv(data_dir + "halfvarson/mapping.txt", sep = "\t")
mapping = mapping.set_index("sample_name")


best_hits = pd.read_csv(data_dir + "halfvarson/embed/best_hits.tsv", header = None, sep = "\t")
best_hits.columns = ["query_id", "hit_id", "query_seq", "hit_seq", "evalue", "bitscore"]
best_hits = best_hits.set_index('query_id')
keep = [i < 1E-29 for i in best_hits['evalue'] ]
best_hits = best_hits.loc[keep, :]
print(best_hits.shape)
best_hits


#Get only those ASVs that have a match in the embedding set

best_hits = best_hits.loc[[i in otu.columns.values for i in best_hits['query_seq']], :]
otu_use = otu.loc[:, best_hits['query_seq'].values]

#Assign id of best match and reorder columns
otu_use.columns = best_hits.loc[:, 'hit_id']

#Put transformation matrix in order to be dotted with the ASV table
qual_vecs = qual_vecs.loc[otu_use.columns.values, :]

map_onesample = mapping.loc[[not i for i in mapping.duplicated('patientnumber')], :]
map_onesample = map_onesample.loc[[i in otu_use.index.values for i in map_onesample.index.values], :]


# In[11]:


otu_use = otu_use.loc[map_onesample.index, :]
otu_use.shape

map_cd_uc = map_onesample.loc[[i == "CD" or i == "UC" for i in map_onesample["diagnosis_full"]], :]
otu_cd_uc = otu_use.loc[map_cd_uc.index.values, :]



import random
importlib.reload(hf)
import warnings
warnings.filterwarnings('ignore')
def crossValPrediction(otu_use, map_onesample, max_depth = 10, n_estimators = 65, weight = 5):
    folds = 5
    kf = KFold(n_splits = folds, shuffle = True, random_state = 110)
    kf.get_n_splits(otu_use)
    
    auc_crossVal = []
    f1_crossVal = []
    feat_imp_crossVal = []
    auc_prec_crossVal = []
    i = 0
    for train_index, val_index in kf.split(otu_use):
        print(val_index)
        #plt.subplot(2,2, i + 1)
        otu_train = otu_use.iloc[train_index, :]
        otu_val = otu_use.iloc[val_index, :]
        map_train = map_onesample.iloc[train_index, :]
        map_val = map_onesample.iloc[val_index, :]
        y_train = map_train['diagnosis_full'].values == "UC"
        y_val = map_val['diagnosis_full'].values == "UC"
        
        auc, auc_train, fpr, tpr, avg_prec, f1, feat_imp = hf.predictIBD(otu_train, y_train, otu_val, y_val,
                  max_depth = max_depth, n_estimators = n_estimators, weight = weight, plot = True, plot_pr = True, feat_imp = True)
        
        auc_crossVal.append(auc)
        auc_prec_crossVal.append(avg_prec)
        f1_crossVal.append(f1)
        feat_imp_crossVal.append(feat_imp)
        
        i = i + 1
    return(auc_crossVal,auc_prec_crossVal, f1_crossVal, feat_imp_crossVal)

# In[61]:




f = plt.figure(figsize=(15,5))
auc_crossVal, auc_prec_crossVal, f1_crossVal, feat_imp_crossVal = crossValPrediction(hf.asinh(otu_cd_uc), map_cd_uc, max_depth = 10, n_estimators = 65, weight = 5)
print("Asinh")
print(np.mean(auc_crossVal))
print(np.mean(auc_prec_crossVal))
f.savefig(fig_dir + "curves_halfvarson_asin.pdf")


embedded = pd.DataFrame(np.dot(hf.asinh(otu_cd_uc), qual_vecs))
embedded.columns = qual_vecs.columns.values

f = plt.figure(figsize=(15,5))
auc_crossVal, auc_prec_crossVal, f1_crossVal, feat_imp_crossVal = crossValPrediction(embedded, map_cd_uc, max_depth = 10, n_estimators = 170,  weight = 5)
print("Embed")
print(np.mean(auc_crossVal))
print(np.mean(auc_prec_crossVal))
f.savefig(fig_dir + "curves_halfvarson_embed.pdf")



f = plt.figure(figsize=(15,5))
from sklearn.decomposition import PCA
pca = PCA(n_components= qual_vecs.shape[1])
pca.fit(hf.asinh(otu_cd_uc))
otu_pca = pd.DataFrame(pca.transform(hf.asinh(otu_cd_uc)))
auc_crossVal, auc_prec_crossVal, f1_crossVal, feat_imp_crossVal = crossValPrediction(otu_pca, map_cd_uc,
                                                 max_depth = 10, n_estimators = 110, weight = 5)
print("PCA")
print(np.mean(auc_crossVal))
print(np.mean(auc_prec_crossVal))
f.savefig(fig_dir + "curves_halfvarson_pca.pdf")
