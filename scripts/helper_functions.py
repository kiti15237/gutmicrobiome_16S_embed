import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import confusion_matrix
import math
import copy
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.model_selection import train_test_split


from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc

def asinh(otu):
    return(np.arcsinh(otu))

def log_normalize(otu):
    return(np.log(otu + 1))

def is_number(s):
    try:
        if np.isnan(float(s)):
            return False
        else: 
            return True
    except ValueError:
        return False

def getDataAG(otu_file, mapping_file, qual_vec_file, test_samples_file, number_criteria, cat_criteria):
    
    otu = pd.read_csv(otu_file, sep = "\t", index_col= 0)
    otu.index = otu.index.map(str)
    #Delete samples that have no taxa
    #otu = otu.loc[:, otu.sum(axis=0) != 0]
    otu.head()
    print("Original data dimensions")
    print("Taxa: " + str(otu.shape[0]) + "  Samples: " + str(otu.shape[1]))


    
    #Format quality vector matrices
    qual_vecs = pd.read_csv(qual_vec_file, sep = " ", index_col = 0, header=None, dtype = {0:str})
    otu_clean = otu.loc[[i in qual_vecs.index.values for i in otu.index.values]] #Keep taxa if present in quality vectors
    otu_clean = otu_clean.reindex(sorted(otu_clean.index.values), axis = 0) #Sort taxa
    qual_vecs_clean = qual_vecs.loc[[i in otu_clean.index.values for i in qual_vecs.index.values]] #keep taxa if present in otu
    qual_vecs_clean = qual_vecs_clean.reindex(sorted(qual_vecs_clean.index.values), axis = 0) #Sort taxa
    print("Filter for taxa present in embeddings")
    print("Taxa: " + str(otu_clean.shape[0]) + "  Samples: " + str(otu_clean.shape[1]))
    

    mapping = pd.read_csv(mapping_file, sep = "\t", index_col=0)
    map_clean = mapping.loc[otu_clean.columns.values] #Keep samples if present in otu
    map_clean = map_clean.reindex(sorted(map_clean.index.values))
    otu_clean = otu_clean.reindex(sorted(otu_clean.columns.values), axis = 1)

    
    
    keep = [True] * map_clean.shape[0]
    print("Samples originally: " + str(sum(keep)))
    for criteria in cat_criteria:
        #Keep sample if it has desired metadata available
        keep_tmp = [( (i != "Unknown") and (i != "Unspecified") and (i!="other" ) and (i != "unspecified") and (isinstance(i, str)) )  for i in map_clean[criteria]] 
   
        keep = [(i and j) for (i,j) in zip(keep, keep_tmp)]
    print("Samples after categorical filter: " + str(sum(keep)))

    for criteria in number_criteria:
        keep_tmp = [is_number(i) for i in map_clean[criteria]]
        keep = [i and j for (i,j) in zip(keep, keep_tmp)] 
    print("Samples after numerical filter: " + str(sum(keep)))

        
    otu_keep = otu_clean.loc[:, keep]
    map_keep = map_clean.loc[keep, cat_criteria + number_criteria]   
    otu_keep = otu_keep.T

    print("Filter for desired metadata present")
    print("Samples: " + str(otu_keep.shape[0]) + "  Taxa: " + str(otu_keep.shape[1]))

    #Make train/test set
    f = open(test_samples_file, "rb")
    test_samples = pickle.load(f)
    f.close

    test_samples = test_samples[[i in otu_keep.index.values for i in test_samples]] #Only include ids that we haven't dropped in previous steps
    otu_train = otu_keep.loc[[not(i in test_samples) for i in otu_keep.index.values], :]
    otu_test = otu_keep.loc[[i in test_samples for i in otu_keep.index.values], :]

    map_train = map_keep.loc[[not(i in test_samples) for i in map_keep.index.values], :]
    map_test = map_keep.loc[[i in test_samples for i in map_keep.index.values], :]

    print(otu_train.shape)
    print(map_train.shape)
    print(otu_test.shape)
    print(map_test.shape)


    return(otu_train, otu_test, qual_vecs_clean, map_train, map_test)


def getDataFreshwater(qual_vec_type = "500"):
    otu = pd.read_csv("data/silva/freshwater/otu_filtered_freshwater.csv", sep = "\t", index_col=0)
    print(otu.shape)
    otu.head()
    #Delete samples that have no taxa
    otu = otu.loc[:, otu.sum(axis=0) != 0]
    
    #Read quality vector data
    if qual_vec_type == "100":
        qual_vecs = pd.read_csv("embeddings/silva/glove_emb_freshwater_100.txt", sep = " ", index_col = 0, header=None)
    if qual_vec_type == "500":
        qual_vecs = pd.read_csv("embeddings/silva/glove_emb_freshwater_2perc_500.txt", sep = " ", index_col = 0, header=None)
    print(qual_vecs.shape)
    qual_vecs.head()
    
    #Match taxa present in quality vectors with those in otu table
    bools_drop = [i not in qual_vecs.index.values for i in otu.index.values]
    drop = otu.index.values[bools_drop]
    otu_drop = otu.drop(drop, axis = 0)
    print("OTU DROP SHAPE: " + str(otu_drop.shape))
    
    qual_vecs = qual_vecs.drop("<unk>", axis = 0)
    qual_vecs_sort = qual_vecs.reindex(sorted(qual_vecs.index.values), axis = 0)
    otu = otu_drop.reindex(sorted(otu_drop.index.values), axis = 0) #Organize otu rows to match taxa in qual_vecs

    ntaxa = qual_vecs.shape[0]
    print(ntaxa)
    bools_correct = [qual_vecs_sort.index.values[i] == otu.index.values[i] for i in range(ntaxa)]

    if (np.sum(bools_correct) == qual_vecs_sort.shape[0]) and (np.sum(bools_correct) == otu.shape[0]):
        print("Safe to continue")
    else:
        print("STOP! Something is wrong.")
        
        
    #Read mapping data
    mapping = pd.read_csv("data/emp_qiime_mapping_release1.tsv.csv", sep = ",", index_col=0, encoding='latin1')
    print("Mapping has shape: " + str(mapping.shape))

    map_filt = mapping.loc[mapping['empo_3'] == "Water (non-saline)"]
    print("After selecting for biome, mapping has shape " + str(map_filt.shape))

    bools = [i in otu.columns.values for i in map_filt.index.values]
    map_filt = map_filt.loc[bools]
    print("After selecting just the samples present in the otu table: " + str(map_filt.shape))

    otu_sort = otu.reindex(sorted(otu.columns.values), axis = 1)
    map_filt_sort = map_filt.reindex(sorted(map_filt.index.values), axis = 0)
    nsamples = map_filt.shape[0]
    bools_correct = [map_filt_sort.index.values[i] == otu_sort.columns.values[i] for i in range(nsamples)]
    print("After rearranging, we have " + str(np.sum(bools_correct)) + " matching samples")
    if (np.sum(bools_correct) == map_filt_sort.shape[0]) and (np.sum(bools_correct) == otu_sort.shape[1]):
        print("Safe to continue")
    else:
        print("STOP! Something is wrong.")
        
        
    #temperature
    #phosphate
    #ammonia
    #map_filt_sort.loc[map_filt_sort['temperature_deg_c']]
    bools =  ~map_filt_sort['temperature_deg_c'].isnull() 
    map_temp = map_filt_sort.loc[~map_filt_sort['temperature_deg_c'].isnull()]
    otu_temp = otu_sort.loc[:, bools]
    bools_correct = [map_temp.index.values[i] == otu_temp.columns.values[i] for i in range(map_temp.shape[0])]
    print("We will be working with " + str(np.sum(bools_correct)) + " samples that have temperature information")
    
    
    #One final check after all transformations
    qual_vecs = qual_vecs_sort
    bools_correct = [qual_vecs.index.values[i] == otu_temp.index.values[i] for i in range(ntaxa)]
    if (np.sum(bools_correct) == qual_vecs_sort.shape[0]) and (np.sum(bools_correct) == otu.shape[0]):
        print("Safe to continue")
    else:
        print("STOP! Something is wrong.")
        
    pd_qual_vecs = pd.DataFrame(qual_vecs)
    otu = otu_temp.T
    
    otu_train = otu.loc[map_temp.study_id != 1883, :]
    otu_test = otu.loc[map_temp.study_id == 1883, :]
    map_train = map_temp.loc[map_temp.study_id != 1883, :]
    map_test = map_temp.loc[map_temp.study_id == 1883, :]
    #map_temp = map_temp.loc[map_temp.study_id != 1041, :]
    return(otu_train, otu_test, pd_qual_vecs, map_train, map_temp)


def normalize(otu):
    #Normalize
    sample_sums = otu.sum(axis=1)
    otu_norm = otu.div(sample_sums, axis=0)
    return(otu_norm)

def biofilter(otu):
    #Filter for useful taxa
    file = open('feature_selection/taxa_lowphy_highcos.obj', 'rb')
    taxa_lowphy_highcos = pickle.load(file)
    file.close()
    otu_use = otu[list(taxa_lowphy_highcos)]
    return(otu_use)

def embed(otu, qual_vecs):
    qual_vecs_use = qual_vecs.loc[list(otu.columns.values)]
    df = pd.DataFrame(np.dot(otu, qual_vecs_use), index = otu.index.values)
    return(df)

def test():
    print("Hello world")
    
    
def plotLineOfBestFit(xi, y, title = ""):
    #plt.figure(figsize=(15,5))
    print(np.max(y))

    #y = y.values
    slope, intercept, r_value, p_value, std_err = stats.linregress(xi,y)
    line = slope*xi+intercept
    plt.plot(xi,y,'o', color = "#AC76C8")
    plt.plot(xi, line, color = "orange")
    
    #perfect line
    #plt.plot(
    
    #plt.xlabel("Predicted Value")
    #plt.ylabel("True Value")
    #plt.title(title)
    plt.ylim((np.min(y),np.max(y)))
    plt.xlim((np.min(y), np.max(y)))
    print("Slope: " + str(slope))
    print("R value: " + str(r_value))
    
    
def testModel(X_train, y_train, X_test, y_test, model, title):
    if model == "linreg":
        m = linear_model.LinearRegression()
    if model == "svm":
        m = svm.SVR(gamma = "scale")
    if model == "rf":
        if is_number(y_train[0]):
            m = RandomForestRegressor(max_depth= 12, random_state=0, n_estimators=140)
        else:
            m = RandomForestClassifier(max_depth= 12, random_state=0, n_estimators=140)
    m.fit(X_train, y_train)
    preds = m.predict(X_test)
    
    plotLineOfBestFit(preds, y_test, title)

    error = np.mean([np.absolute(y_test[i] - preds[i]) for i in range(len(preds))])
    mse = np.mean([np.square(y_test[i] - preds[i]) for i in range(len(preds))])
    print("Linear Error" + str(error) + "   MSE: " + str(mse))
    return(m, error, mse)
    #plotLineOfBestFit(y_test, preds, title = title)
#testModel(X_train, y_train, X_test, y_test, model = "linreg", title = "")



from sklearn.utils.multiclass import unique_labels
def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=True,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax


def testModelClass(X_train, y_train, X_test, y_test, model = 'rf', weights={0:1, 1:1}, title = ""):
    m = OneVsRestClassifier(RandomForestClassifier(max_depth= 12, random_state=0, n_estimators=140, class_weight = weights))
    m.fit(X_train, y_train)
    preds = m.predict(X_test)
    plot_confusion_matrix(y_test, preds, classes = np.unique(y_test), title = title)
    
    
    
def makeNumeric(var, dictionary, map_train, map_test, map_train_correct, map_test_correct):
    map_train_correct[var] = [dictionary[i] for i in map_train[var]]
    map_test_correct[var] = [dictionary[i] for i in map_test[var]]
    return(map_train_correct, map_test_correct)

def makeMappingNumeric(map_train, map_test, number_criteria, cat_criteria):
    map_train_correct = copy.deepcopy(map_train)
    map_test_correct = copy.deepcopy(map_test)
    
    for num_var in number_criteria:
        print(num_var)
        map_train_correct[num_var] = [float(i) for i in map_train[num_var]]
        map_test_correct[num_var] = [float(i) for i in map_test[num_var]]
  
    
    #Deal with variables that have to do frequency
    freq_dict = {"Never": 0, "Rarely (a few times/month)":1, "Occasionally (1-2 times/week)": 2,
             "Regularly (3-5 times/week)":3, "Daily":4, "Rarely (less than once/week)" : 1}
    
    freq_vars = ["EXERCISE_FREQUENCY", "ONE_LITER_OF_WATER_A_DAY_FREQUENCY", "SEAFOOD_FREQUENCY", "PROBIOTIC_FREQUENCY", 
             "OLIVE_OIL", "FRUIT_FREQUENCY", "SUGAR_SWEETENED_DRINK_FREQUENCY", "MILK_CHEESE_FREQUENCY", "RED_MEAT_FREQUENCY",
            "MEAT_EGGS_FREQUENCY", "VEGETABLE_FREQUENCY"]
    keep = [(i in cat_criteria) for i in freq_vars]
    freq_vars = np.array(freq_vars)[keep]
    for var in freq_vars:
        print(var)
        map_train_correct, map_test_correct = makeNumeric(var, freq_dict, map_train, map_test, map_train_correct, map_test_correct)

    #Deal with sleep
    if "SLEEP_DURATION" in cat_criteria:
        print("SLEEP_DURATION")
        sleep_dict = {"Less than 5 hours":1, "5-6 hours":2, "6-7 hours":3, "7-8 hours":4, "8 or more hours":5 }
        map_train_correct, map_test_correct = makeNumeric("SLEEP_DURATION", sleep_dict, map_train, map_test, map_train_correct, map_test_correct)

    if "SEX" in cat_criteria:
        print("SEX")
        sex_dict = {"male": 0, "female":1}
        map_train_correct, map_test_correct = makeNumeric("SEX", sex_dict, map_train, map_test, map_train_correct, map_test_correct)

    if "IBD" in cat_criteria:
        print("IBD")
        ibd_dict = {'I do not have this condition':0, 'Self-diagnosed':1, 
                    'Diagnosed by a medical professional (doctor, physician assistant)':1,
                    'Diagnosed by an alternative medicine practitioner':1}
        map_train_correct, map_test_correct = makeNumeric("IBD", ibd_dict, map_train, map_test, map_train_correct, map_test_correct)
    
    if "GLUTEN" in cat_criteria:
        print("GLUTEN")
        gluten_dict = {"No":0, 'I do not eat gluten because it makes me feel bad':1, 'I was diagnosed with celiac disease':1,
                       'I was diagnosed with gluten allergy (anti-gluten IgG), but not celiac disease':1}
        map_train_correct, map_test_correct = makeNumeric("GLUTEN", gluten_dict, map_train, map_test, map_train_correct, map_test_correct)
    return(map_train_correct, map_test_correct)

def predictIBD(X_train, y_train, X_test, y_test, graphTitle = "", max_depth = 12, n_estimators = 140, plot = True):
    weight_list = [1, 5, 20, 100]
    auc_list = []
    auc_train_list = []
    for weight in weight_list:
        weights = {0:1, 1:weight}
        m = OneVsRestClassifier(RandomForestClassifier(max_depth= max_depth, random_state=0, n_estimators= n_estimators, class_weight = weights))
        m.fit(X_train, y_train)
        probs = m.predict_proba(X_test)
        probs_train = m.predict_proba(X_train)

        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y_test, probs[:, 1])
        fpr_train, tpr_train, thresholds_train = roc_curve(y_train, probs_train[:, 1])
                                                     
        roc_auc = auc(fpr, tpr)
        roc_auc_train = auc(fpr_train, tpr_train)
        auc_list.append(roc_auc)
        auc_train_list.append(roc_auc_train)
        if plot:
            plt.plot(fpr, tpr, lw=2, alpha=0.3,label='Weight %d (AUC = %0.2f)' % (weight, roc_auc))
    if plot:
        plt.legend(loc="lower right")
        x = np.linspace(0, 1, 10)
        plt.plot(x, x)
        plt.title(graphTitle)
    return(auc_list, auc_train_list)

def getPCAReduced(X_train, X_val, X_test, components = 500):
    pca = PCA(n_components= components)
    pca.fit(X_train)
    X_train_pca = pca.transform(X_train)
    X_val_pca = pca.transform(X_val)
    X_test_pca = pca.transform(X_test)
    return(X_train_pca, X_val_pca, X_test_pca)

def plotPCA(table, otu_raw, components):
    pca = PCA(n_components= components)
    pca = pca.fit(table)
    table_pca = pca.transform(table)
    table_pca = table_pca / np.max(table_pca)
    df = pd.DataFrame(table_pca, index = table.index.values)
    sample_sums_table = otu_raw.sum(axis = 1)
    plt.scatter(df.iloc[:,0], df.iloc[:,1], c = sample_sums_table, cmap='viridis')
    plt.colorbar()
    plt.xlabel(pca.explained_variance_ratio_[0])
    plt.ylabel(pca.explained_variance_ratio_[1])
    
    
    

