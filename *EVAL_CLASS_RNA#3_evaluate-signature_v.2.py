### import modules
import urllib
import numpy as np
import pylab as plt
import pandas as pd
import multiprocessing as mp
import subprocess
import sys, os, inspect
from scipy import stats
sys.path.append('/Users/dressler/Documents/Core/_PythonModules')
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
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

#########################
data_folder = "*PAULA"
top = int(raw_input("... select top X features: "))
classes = int(raw_input("... classes: "))
#########################
to_predict = pd.read_excel(data_folder + "/CLASSIFIER/RNA/RES_RNA_binary.xlsx", index_col=0).transpose()
data = pd.read_excel(data_folder + "/CLASSIFIER/RNA/RES_binary.xlsx", index_col=0).transpose() #contains clusters, IDs removed as they are indices
candidates = pd.read_excel(data_folder + "/CLASSIFIER/RNA/RES_features-RFECV_frequencies.xlsx")
candidates = candidates.sort_values(by='frequency', ascending=False)
candidates = candidates.reset_index(drop=True)
candidates = candidates["feature"][:top]
print("... selected candidates: %s"%candidates)


y = data.iloc[:,-1]
if classes == 3:
    y = y.replace([1.,2.,3.,4.],[1.,1.,1.,2.])
    y = y.replace([0.,1.,2.],[0,1,2])
#label_encoder = LabelEncoder()
#y = label_encoder.fit_transform(y).astype('int32')
X = data.iloc[:,:-1]
X = X[candidates]
to_predict = to_predict[candidates]



vs.style('science')
colors = vs.colors_indexed
fig = plt.figure("Classifier", [8,8])
ax = fig.add_subplot(221)
ax.set_ylim(0,1.)
ax.set_xticklabels([])
ax.set_ylabel("Accuracy [1]")
def plot(ax,x,y,color,len_all):
    """Plot boxplot with swarmplot"""
    vp = ax.violinplot([y], positions=[x], showmeans=False, showmedians=False, showextrema=False, widths = .9)
    for param in vp['bodies']:
        param.set_alpha(0.25)
        param.set_facecolor(color)
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
    sns.swarmplot(np.array([x] * len(y)), y, order=range(len_all), color=color, alpha=.5, ax=ax, size=3)
    ax.set_xticks(range(len_all))
    ax.set_xlim(-1,len_all)
    ax.grid(False)
    pass

accuracies = []
errors = []
features = X.columns
runs = 100
for k in range(runs):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2)
    

    svc=LinearSVC() # The default kernel used by SVC is the gaussian kernel kernel="linear"
    #print X_train, y_train
    svc.fit(X_train, y_train)
    prediction = svc.predict(X_test)

    #prediction = f.predict(X_test)
    #print("... selected features: %s"%features_select)
    #print prediction
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
print features
plot(ax, 0., accuracies, colors[0], 2)
plot(ax, 1., 1-np.array(errors), colors[0], 2)
ax.set_xticklabels(["accuracy","accuracy I/III"], ha='right', rotation=45.)
res = pd.DataFrame({"accuracy":accuracies, "relevant error":errors})
res.to_excel(data_folder + "/CLASSIFIER/RNA/RES_features-RFECV_evaluate.xlsx", index=False)
fig.savefig(data_folder + "/CLASSIFIER/RNA/PLOT_features-RFECV_avaluate.png", dpi=300)

### predict
#vc=SVC()
#svc.fit(X, y)
svc=SVC(probability=True)
svc.fit(X, y)
prediction = svc.predict_proba(to_predict)
print prediction

#prediction = pd.DataFrame({"sample":to_predict.index, "cluster":prediction})
prediction = pd.DataFrame(prediction)
print len(prediction.idxmax(axis=1).values)
prediction["cluster"] = prediction.idxmax(axis=1).values
prediction["cluster_prob"] = prediction[[x for x in prediction.columns if x != "cluster"]].max(axis=1).values
prediction["sample"] = to_predict.index

#prediction["sample"] = prediction["sample"].str.split(",").str[0]
#prediction["sample"] = prediction["sample"].str.replace("-", "")
prediction.to_excel(data_folder + "/CLASSIFIER/RNA/RES_prediction.xlsx", index=False)
sys.exit()
prediction = svc.predict(to_predict)
print prediction

pd.DataFrame({"sample":to_predict.index, "cluster":prediction}).to_excel(data_folder + "/CLASSIFIER/RNA/RES_prediction.xlsx", index=False)


plt.show()