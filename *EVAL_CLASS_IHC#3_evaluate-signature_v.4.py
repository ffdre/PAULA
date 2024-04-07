### import modules
import urllib
import numpy as np
import pylab as plt
import pandas as pd
import multiprocessing as mp
import subprocess
import sys, os, inspect
from scipy import stats
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
sys.path.append(file_directory + '/_modules')
from statsmodels.stats import multitest
import math
import visualize as vs
import statistics as st
import string
import random
import subprocess
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
import random
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import confusion_matrix
from sklearn.feature_selection import SelectKBest, chi2,mutual_info_classif
from sklearn.feature_selection import RFECV
from sklearn.svm import SVR
import statsmodels.formula.api as sm
plt.show()
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

#########################
data_folder = "*PAULA"
#top = int(raw_input("... select top X features: "))
top = 2
samples_excluded = []#["U-BEG-81_T","U-CAC-83_T","U-CUA-20_T","U-CUM-84_T","U-PUJ-07_T","U-QAH-42_T", "U-PEC-04_T"]
exclude_features = [
                    "CD99 antigen OS=Homo sapiens (Human) OX=9606 GN=CD99 PE=1 SV=1",
                    "Collagen alpha-1(IV) chain OS=Homo sapiens (Human) OX=9606 GN=COL4A1 PE=1 SV=4",
                    ]
#########################
to_predict = pd.read_excel(data_folder + "/CLASSIFIER/IHC/*RES_IHC-cohort_binary.xlsx", index_col=0)
to_predict = to_predict[[x for x in to_predict.columns if "area" not in x and x not in samples_excluded]]
new_indices = []
for x in to_predict.columns:
    #x = x[:-5]
    if "_T" not in x: new_indices.append(x+"_T")
    else: new_indices.append(x)
to_predict.columns = new_indices
to_predict = to_predict.transpose()
to_predict = to_predict.dropna()
for column in to_predict.columns:
    if "ollag" not in column:
        temp = []
        low = np.percentile(to_predict[column].values,25.)
        up = np.percentile(to_predict[column].values,75.)
        for x in to_predict[column].values:
            if x <= low: temp.append(-1)
            elif x >= up: temp.append(1)
            else: temp.append(0)
        to_predict[column] = temp



to_predict.to_excel(data_folder + "/CLASSIFIER/IHC/CTRL.xlsx")
data = pd.read_excel(data_folder + "/CLASSIFIER/IHC/RES_binary.xlsx", index_col=0).transpose() #contains clusters, IDs removed as they are indices
data_ori = pd.read_excel(data_folder + "/CLASSIFIER/IHC/RES_original.xlsx", index_col=0).transpose() #contains clusters, IDs removed as they are indices
candidates = pd.read_excel(data_folder + "/CLASSIFIER/IHC/RES_features-RFECV_frequencies.xlsx")
#candidates = candidates.sort_values(by='frequency', ascending=False)
candidates = candidates.reset_index(drop=True)
candidates = [x for x in candidates["feature"][:top] if x not in exclude_features]
#candidates = candidates[[x for x in candidates.columns]]
print("... selected candidates: %s"%candidates)

#data.to_excel(data_folder + "/CLASSIFIER/IHC/temp_data.xlsx")
data["Cluster"] = data["Cluster"].replace([1.,2.,3.,4.],[1.,1.,1.,2.])
y = data.iloc[:,-1]

#y = y.replace([0.,1.,2.],[0,1,2])
#y.to_excel(data_folder + "/CLASSIFIER/IHC/temp_y.xlsx")
#label_encoder = LabelEncoder()
#y = label_encoder.fit_transform(y).astype('int32')
X = data_ori.loc[~data_ori.index.isin(to_predict.index.values),:]
X = X[candidates]

for column in X.columns:
    if "ollag" not in column:
        temp = []
        low = np.percentile(X[column].values,25.)
        up = np.percentile(X[column].values,75.)
        for x in X[column].values:
            if x <= low: temp.append(-1)
            elif x >= up: temp.append(1)
            else: temp.append(0)
        X[column] = temp



vs.style('science')
colors = vs.colors_indexed
fig = plt.figure("Classifier", [8,8])
ax = fig.add_subplot(221)
ax.set_ylim(0,1.)
ax.set_xticklabels([])
ax.set_ylabel("Accuracy [1]")
def plot(ax,x,y,color,len_all):
    """Plot boxplot with swarmplot"""
    #vp = ax.violinplot([y], positions=[x], showmeans=False, showmedians=False, showextrema=False, widths = .9)
    #for param in vp['bodies']:
    #    param.set_alpha(0.25)
    #    param.set_facecolor(color)
    ax.boxplot([y], positions=[x], showfliers=False, patch_artist=True, widths=.45,
               boxprops=dict(facecolor=color, alpha=.25, linewidth=.5),
               capprops=dict(color=color, alpha=.5),
               whiskerprops=dict(color=color, alpha=.5),
               flierprops=dict(color="black", markeredgecolor=None, markersize=5., marker="."),
               medianprops=dict(color="grey"))
    if np.median(y) < np.mean(y):
        ax.plot([x-.2, x+.2],[np.mean(y), np.mean(y)], color=color, lw=2., alpha=.5)
    else:
        ax.plot([x-.2, x+.2],[np.mean(y), np.mean(y)], color=color, lw=2., alpha=.5)
    #sns.swarmplot(np.array([x] * len(y)), y, order=range(len_all), color=color, alpha=.5, ax=ax, size=3)
    ax.set_xticks(range(len_all))
    ax.set_xlim(-1,len_all)
    ax.grid(False)
    pass

accuracies = []
errors = []
features = X.columns
runs = 100
ids = X.index.values
res_samples = {}
for id in ids:
    res_samples[id] = []
for k in range(runs):
    X_train, X_test, y_train, y_test, id_train, id_test = train_test_split(X, y[X.index.values], ids, test_size=.2)
    svc=SVC() # The default kernel used by SVC is the gaussian kernel kernel="linear"
    #print X_train, y_train
    svc.fit(X_train, y_train)
    prediction = svc.predict(X_test)
    for sample, pred, true in zip(id_test,prediction,y_test):
        print sample, pred, true
        res_samples[sample].append(np.abs(pred-true))
    #prediction = f.predict(X_test)
    #print("... selected features: %s"%features_select)
    #print y_test, prediction
    cm = confusion_matrix(y_test, prediction)
    print cm
    sum = 0
    sum_margin = cm[0][-1]
    sum_margin += cm[-1][0]
    #print float(X_test.shape[0])
    false_margin = sum_margin / float(X_test.shape[0])
    for i in range(cm.shape[0]):
        sum += cm[i][i]

    accuracy = sum / float(X_test.shape[0])
    print("...accuracy was: %s"%accuracy)
    print("...relevant error was: %s"%false_margin)
    accuracies.append(accuracy)
    errors.append(false_margin)
if 1 == 1:
    print features
    res_samples_eval = pd.DataFrame()
    res_samples_eval["legend"]=["cluster","accuracy", "error neighboring cluster", "error I/III"]
    #print y, y_test, prediction
    for sample, cluster in zip(ids,y):
        temp = res_samples[sample]
        res_samples_eval[sample] = [cluster,temp.count(0)/float(len(temp)), temp.count(1)/float(len(temp)), temp.count(2)/float(len(temp))]
    res_samples_eval.transpose().to_excel(data_folder + "/CLASSIFIER/IHC/RES_performance-across-samples.xlsx")
    plot(ax, 0., accuracies, colors[0], 2)
    plot(ax, 1., 1-np.array(errors), colors[0], 2)
    ax.set_xticklabels(["accuracy","accuracy I/III"], ha='right', rotation=45.)
    res = pd.DataFrame({"accuracy":accuracies, "relevant error":errors})
    res.to_excel(data_folder + "/CLASSIFIER/IHC/RES_features-RFECV_evaluate.xlsx", index=False)
    fig.savefig(data_folder + "/CLASSIFIER/IHC/PLOT_features-RFECV_evaluate.png", dpi=300)
    print("... selected candidates: %s"%candidates)
#plt.show()

### predict
#plt.show()
svc=SVC(probability=True, kernel='linear')
svc.fit(X.loc[~X.index.isin(to_predict.index.values),:], y[~data.index.isin(to_predict.index.values)])
#svc.fit(X, y)
prediction = svc.predict_proba(to_predict)
prediction = pd.DataFrame(prediction)

cm = confusion_matrix(data.loc[to_predict.index.values, "Cluster"].values, prediction.idxmax(axis=1).values)
print prediction
print cm
sum = 0
sum_margin = cm[0][-1]
sum_margin += cm[-1][0]
#print float(X_test.shape[0])
false_margin = sum_margin / float(to_predict.shape[0])
for i in range(cm.shape[0]):
    sum += cm[i][i]

accuracy = sum / float(to_predict.shape[0])
print("...accuracy was: %s"%accuracy)
print("...relevant error was: %s"%false_margin)


#prediction = pd.DataFrame({"sample":to_predict.index, "cluster":prediction})

print len(prediction.idxmax(axis=1).values)

prediction["cluster_pred"] = prediction.idxmax(axis=1).values
#prediction["cluster_true"] = data.loc[to_predict.index.values, "Cluster"].replace([1.,2.,3.,4.],[1.,1.,1.,2.]).values
prediction["cluster_true"] = data.loc[to_predict.index.values, "Cluster"].values
prediction["cluster_prob"] = prediction[[x for x in prediction.columns if x != "cluster_pred" and x != "cluster_true"]].max(axis=1).values

prediction["sample"] = to_predict.index

correct = len(prediction["cluster_true"].where(prediction["cluster_true"] == prediction["cluster_pred"]).dropna())
all = len(prediction["cluster_true"])

print("... correct: %s\t... all: %s\t... share: %s %%"%(correct, all, correct/float(all)*100))

#prediction["sample"] = prediction["sample"].str.split(",").str[0]
#prediction["sample"] = prediction["sample"].str.replace("-", "")
prediction.to_excel(data_folder + "/CLASSIFIER/IHC/*RES_cohort-prediction.xlsx", index=False)


sys.exit()
svc=SVC()
svc.fit(X, y)
prediction = svc.predict(to_predict)
print prediction

pd.DataFrame({"sample":to_predict.index, "cluster":prediction}).to_excel(data_folder + "/CLASSIFIER/IHC/*RES_cohort-prediction.xlsx", index=False)


plt.show()