
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import helper_functions as hf
import pickle
import matplotlib.pyplot as plt
import importlib
import math
import copy
import re
from sklearn import preprocessing
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
import random
importlib.reload(hf)
import warnings
warnings.filterwarnings('ignore')



data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"
fig_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/figures/"

qual_vecs, embed_ids, embed_seqs = hf.getQualVecs(data_dir)



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



otu_use = otu_use.loc[map_onesample.index, :]
otu_use.shape

map_cd_uc = map_onesample.loc[[i == "CD" or i == "UC" for i in map_onesample["diagnosis_full"]], :]
otu_cd_uc = otu_use.loc[map_cd_uc.index.values, :]
y = map_cd_uc['diagnosis_full'].values == "UC"


f = plt.figure(figsize=(15,5))
auc_crossVal, auc_prec_crossVal, f1_crossVal, feat_imp_asin = hf.crossValPrediction(hf.asinh(otu_cd_uc), y, max_depth = 2, n_estimators = 50, weight = 20, folds = 5)
print("Asinh")
print(otu_cd_uc.shape)
print(np.mean(auc_crossVal))
print(np.mean(auc_prec_crossVal))
f.savefig(fig_dir + "curves_halfvarson_asin_tmp.pdf")


embedded = pd.DataFrame(np.dot(hf.asinh(otu_cd_uc), qual_vecs))
embedded.columns = qual_vecs.columns.values

f = plt.figure(figsize=(15,5))
auc_crossVal, auc_prec_crossVal, f1_crossVal, feat_imp_embed = hf.crossValPrediction(embedded, y, max_depth = 2, n_estimators = 50,  weight = 20, folds = 5)
print("Embed")
print(embedded.shape)
print(np.mean(auc_crossVal))
print(np.mean(auc_prec_crossVal))
f.savefig(fig_dir + "curves_halfvarson_embed_tmp.pdf")



f = plt.figure(figsize=(15,5))
from sklearn.decomposition import PCA
pca = PCA(n_components= qual_vecs.shape[1])
pca.fit(hf.asinh(otu_cd_uc))
otu_pca = pd.DataFrame(pca.transform(hf.asinh(otu_cd_uc)))
auc_crossVal, auc_prec_crossVal, f1_crossVal, feat_imp_pca = hf.crossValPrediction(otu_pca, y,
                                                 max_depth = 2, n_estimators = 50, weight = 20, folds = 5)
print("PCA")
print(otu_pca.shape)
print(np.mean(auc_crossVal))
print(np.mean(auc_prec_crossVal))
f.savefig(fig_dir + "curves_halfvarson_pca_tmp.pdf")



####################################################
###############  Feature Importance  ###############
####################################################


feat_imp_df = hf.getFeatImpDf(feat_imp_embed)
pathway_table = pd.read_csv(data_dir + "pathways/property_pathway_dict.txt",
                            sep= "\t", dtype= {"pathway_id": 'object'})
pathway_table = pathway_table.set_index('dim')
tmp = pathway_table.loc[feat_imp_df.index.values, :]
feat_imp_df_paths = pd.merge(feat_imp_df, tmp, left_index = True, right_index = True)


max_depth = 2
n_estimators = 50
weight = 20
weights = {0:1, 1:weight}
m = RandomForestClassifier(max_depth= max_depth, random_state=0, n_estimators=n_estimators, class_weight = weights)
m.fit(embedded, y)


##########################################################################
##### Functions to propogate associations in decision trees  #############
##########################################################################

def getLeafInx(estimator):
    leaf_inx = np.where([i == -1 and j == -1 for i,j in zip(estimator.tree_.children_left, estimator.tree_.children_right)])
    return(leaf_inx[0])

def getLabels(estimator):
    value = estimator.tree_.value
    leaf_inx = getLeafInx(estimator)
    labels = []
    for val_set in value[leaf_inx]:
        val_set = np.squeeze(val_set)
        if val_set[0] > val_set[1]:
            val = "CD"
        else:
            val = "UC"
        labels.append(val)
    return(labels)

def associate_one_level(feature, label_left, label_right):
    associations[feature][label_right] += 1
    associations[feature][label_left] -= 1
    
def left_side_propogate(feature, label_left, label_right):
    if label_left == label_right:
        associations[features[0]][label_left] -= 1 #Doesn't matter which label, they're equal
    else:
        associate_one_level(feature, label_left, label_right)

def right_side_propogate(feature, label_left, label_right):
    if label_left == label_right:
        associations[features[0]][label_left] += 1
    else:
        associate_one_level(feature, label_left, label_right)

def getAssociation(associations, i):
    if associations[i]['UC'] > associations[i]['CD']:
        return('UC')
    if associations[i]['UC'] == associations[i]['CD']:
        return('EQ')
    else: 
        return('CD')
    
def getDiffMag(associations, i):
    return(np.abs(associations[i]['UC'] - associations[i]['CD']))



associations = {}

estimators = m.estimators_
i = 0
for estimator in estimators:
    #print(i)
    features = estimator.tree_.feature + 1 # the forest starts labeling it's features at 0, we start at 1
    features = features[features > 0]
    labels = getLabels(estimator)
    leaf_inx = getLeafInx(estimator)
    op1 = [2, 3, 5, 6] #full set of nodes
    op2 = [2,3,4] #no test node right
    op3 = [1, 3, 4] # no test node left
    for feat in features:
        if not feat in associations:
            associations[feat] = {'UC': 0, 'CD':0}

    if np.array_equal(leaf_inx , op1):
        left_side_propogate(features[1], labels[0], labels[1])
        right_side_propogate(features[2], labels[2], labels[3])

    if np.array_equal(leaf_inx , op2):
        left_side_propogate(features[1], labels[0], labels[1])
        associations[features[0]][labels[2]] += 1

    if np.array_equal(leaf_inx , op3):
        right_side_propogate(features[1], labels[1], labels[2])
        associations[features[0]][labels[0]] -=1
        
    i += 1
print(associations)

topics = [int(re.sub("property_100_", "", i)) for i in feat_imp_df_paths.index.values]
association = np.array(['NA'] * len(topics))
diffMag = np.zeros(len(topics))
i = 0
for topic in topics:
    if topic in associations.keys():
        association[i] = getAssociation(associations, topic)
        diffMag[i] = getDiffMag(associations, topic)
    else:
        association[i] = "NA"
    i +=1
    
feat_imp_df_paths.insert(3, "Association", association)
feat_imp_df_paths.insert(4, "Diff Num. Trees Associated", diffMag)
feat_imp_df_paths
feat_imp_df_paths.to_csv(data_dir + "halfvarson/metabolic_pathways_importance.csv")