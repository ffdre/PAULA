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
from scipy.cluster.hierarchy import dendrogram, linkage
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

#########################
data_folder = "*PAULA"
filter_RNA = False
filter_corr = False
reorder_plot = True
sig = .05
min_corr = .75
clusters = 5
num_features = 1
reverse_log2 = False
zscore = False #must be False if quantile True
quantiles = True
fill_na = False
len_labels = 10
fontsize = 10
fontsize_heatmap = 7
legend = []
features = []
#########################

try: os.mkdir(data_folder + "/CLASSIFIER")
except: pass
try: os.mkdir(data_folder + "/CLASSIFIER/RNA")
except: pass

### load data


data = pd.read_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICSmed-log2.xlsx") #_COMPLETE-DATA
data["Accession"] = pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins.xlsx").loc[:,"Accession"]
data = pd.merge(data, \
    pd.read_excel(data_folder + "/*RES_data_functional.xlsx")[["Accession","Gene Symbol"]], \
    how='left', on=["Accession"])
data["Accession"] = data["Gene Symbol"]
data = data.drop(labels=["Gene Symbol"],axis=1)

if filter_RNA == True:
    RNA = pd.read_excel(data_folder + "/CLASSIFIER/RNA/TCGA_data_RNA_Seq_v2_expression_median.xlsx").iloc[:,:]
    RNA = pd.merge(RNA, data.loc[:,"Accession"], how='inner', left_on=["Hugo_Symbol"], right_on=["Accession"]) #"/PAIR/*RES_data_paired.xlsx"
    RNA = RNA.drop(labels=["Entrez_Gene_Id","Accession"], axis=1)
    RNA.to_excel(data_folder + "/CLASSIFIER/RNA/RNA_data_filtered.xlsx", index=False)

if filter_corr == True:
    RNA = pd.read_excel(data_folder + "/CLASSIFIER/RNA/RNA_data_filtered.xlsx")
    data_corr = pd.read_excel(data_folder + "/CLASSIFIER/RNA/data_correlation.xlsx")
    for i,row in RNA.iterrows():
        try:
            corr = data_corr.where(data_corr["Gene name"] == row["Hugo_Symbol"]).dropna(thresh=1)
            #print corr
            print i
            p = corr["p-value"].values[0]
            r = corr["Spearman correlation"].values[0]
            print p, r
            if p > sig or r < min_corr:
                RNA = RNA.drop(labels=[i], axis=0)
        except:
            RNA = RNA.drop(labels=[i], axis=0)
    RNA.to_excel(data_folder + "/CLASSIFIER/RNA/RNA_data_filtered_positively-correlated.xlsx", index=False)
else:
    RNA = pd.read_excel(data_folder + "/CLASSIFIER/RNA/RNA_data_filtered_positively-correlated.xlsx")

RNA["Accession"] = RNA["Hugo_Symbol"]
RNA = RNA.drop(labels=["Hugo_Symbol"],axis=1)
RNA = RNA.reset_index(drop=True)

data = pd.merge(data, RNA["Accession"], how="right", on=["Accession"])

data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")
print len(data_clinical.iloc[:,0])
data_clinical = data_clinical.dropna(subset = ["cluster k=5"])
print len(data_clinical.iloc[:,0])
data_clinical = data_clinical.where(data_clinical["ID_set"].isin(data.columns)).dropna(thresh=1)
data_clinical = data_clinical.sort_values(by="cluster k=5", ascending=True)





## reverse log2
if reverse_log2 == True:
    data.loc[:,sample_columns] = 2. ** data.loc[:,sample_columns]

## fill NA, prepare sorting, transform to zscore
sorted_samples = RNA.columns[:-1]
sample_columns = ~RNA.columns.isin(['Accession'])

minimum = np.min(np.min(RNA.where(RNA > 0)))
print minimum
RNA.loc[:,sample_columns] = RNA.loc[:,sample_columns].replace(0,minimum)
RNA.loc[:,sample_columns] = np.log2(RNA.loc[:,sample_columns])
for i, row in RNA.iterrows():
    mean = np.nanmean(RNA.iloc[i, sample_columns])
    SD = np.nanstd(RNA.iloc[i, sample_columns], ddof=0)
    if zscore == True:
        RNA.iloc[i, sample_columns] = (RNA.iloc[i, sample_columns] - mean)/float(SD)
        if fill_na == True:
            RNA.iloc[i, sample_columns] = RNA.iloc[i, sample_columns].fillna(0)
    elif quantiles == True:
        up = np.percentile(RNA.iloc[i, sample_columns],75.)
        down = np.percentile(RNA.iloc[i, sample_columns],25.)
        for j, sample in enumerate(sorted_samples):
            if row[sample] >= up: RNA.loc[i, sample] = 1
            elif row[sample] <= down: RNA.loc[i, sample] = -1
            else: RNA.loc[i, sample] = 0
        if fill_na == True:
            RNA.iloc[i, sample_columns] = RNA.iloc[i, sample_columns].fillna(0)
    else:
        if fill_na == True:
            RNA.iloc[i, sample_columns] = RNA.iloc[i, sample_columns].fillna(mean)

RNA[["Accession"]+list(sorted_samples)].to_excel(data_folder + "/CLASSIFIER/RNA/RES_RNA_binary.xlsx", index=False)

sorted_samples = list(data_clinical["ID_set"])
data = data[["Accession"] + sorted_samples]
sample_columns = ~data.columns.isin(['Accession'])
for i, row in data.iterrows():
    mean = np.nanmean(data.iloc[i, sample_columns])
    SD = np.nanstd(data.iloc[i, sample_columns], ddof=0)
    if zscore == True:
        data.iloc[i, sample_columns] = (data.iloc[i, sample_columns] - mean)/float(SD)
        if fill_na == True:
            data.iloc[i, sample_columns] = data.iloc[i, sample_columns].fillna(0)
    elif quantiles == True:
        up = np.percentile(data.iloc[i, sample_columns],75.)
        down = np.percentile(data.iloc[i, sample_columns],25.)
        for j, sample in enumerate(sorted_samples):
            if row[sample] >= up: data.loc[i, sample] = 1
            elif row[sample] <= down: data.loc[i, sample] = -1
            else: data.loc[i, sample] = 0
        if fill_na == True:
            data.iloc[i, sample_columns] = data.iloc[i, sample_columns].fillna(0)
    else:
        if fill_na == True:
            data.iloc[i, sample_columns] = data.iloc[i, sample_columns].fillna(mean)

data_vis = data[sorted_samples]
if reorder_plot == True:
    Y = linkage(data_vis)
    reorder = dendrogram(Y)['leaves']
    data_vis = data_vis.iloc[reorder,:]
data[["Accession"]+sorted_samples].append(pd.DataFrame([["Cluster"]+list(data_clinical["cluster k=5"])], \
                                                       columns=data.columns)).to_excel(data_folder + "/CLASSIFIER/RNA/RES_binary.xlsx", index=False)

### PLOT
fig = plt.figure("expression heatmap",[13.9, 12])
fig.tight_layout()
fig.subplots_adjust(left=.1, right=.9)
GS = gridspec.GridSpec(38, 2, figure=fig, wspace=.0, hspace=0.2, #0.018
                              width_ratios=[0.4, 1])
ax = fig.add_subplot(GS[num_features:,1])

vmax = np.max(np.max(data_vis))
vmin = -1*vmax

cmap = mpl.cm.RdBu_r
cmap.set_bad('grey',1.)
ax.imshow(data_vis, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
ax.grid(False)
ax.set_yticks(range(len(data["Accession"])))
labels = []
for label in data["Accession"]:
    name = label.split(";")[0]
    if len(name) > len_labels:
        name = name[:len_labels] + "..."
    labels.append(name)
ax.set_yticklabels(labels, fontsize=fontsize_heatmap)
features.append(ax)

# colorbar heatmap
cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax = fig.add_subplot(4,120,363)
cb = fig.colorbar(sm,cax=ax, orientation='vertical') #label='fold change [log2]'
cb.ax.tick_params(labelsize=8)
cb.outline.set_linewidth(.0)
ax.set_ylabel("z-score [1]", fontsize=fontsize)

## Features
ax_cb = fig.add_subplot(GS[:,0])
ax_cb.axis('off')

def add_feature(AX, data, name, cmap, names):
    global legend, features
    elements = np.sort(data.unique())
    legend.append((None, name))
    for n, e in zip(names, elements):
        legend.append((mpl.cm.get_cmap(cmap)((e - np.min(elements)) / float(np.max(elements) - np.min(elements))), n))
    legend.append((None, ""))
    AX.imshow([data], cmap=cmap, aspect="auto")
    AX.set_yticks([0])
    AX.set_yticklabels([name], fontsize=fontsize)
    features.append(AX)
    pass

def add_feature_binary(AX, data, name, colors, names):
    global legend, features
    elements = np.sort(data.unique())
    legend.append((None, name))
    for n, e in zip(names, elements):
        legend.append((colors[int(e)], n))
    legend.append((None, ""))
    AX.imshow([data], cmap=mpl.colors.ListedColormap(colors), aspect="auto")
    AX.set_yticks([0])
    AX.set_yticklabels([name], fontsize=fontsize)
    features.append(AX)
    pass

### FEATURES ###

ax = fig.add_subplot(GS[0, 1])
data = data_clinical["cluster k=5"]
name = "Cluster"
names = ["PAULA I", "PAULA II", "PAULA III", "PAULA IV", "PAULA V" ]
cmap = vs.colors_cluster #"Spectral" #"Spectral" Set1
#add_feature(ax, data, name, cmap, names)
#cmap = colors#"Spectral"
add_feature_binary(ax, data, name, cmap, names)


### rest ###
for ax in features:
    ax.grid(False)
    ax.set_frame_on(False)
    ax.tick_params(left=False, right=False)
    ax.set_xticks([])
    ax.set_xticklabels([])


ax_cb.set_ylim(-50,.4)
ax_cb.set_xlim(0,3)
for j,L in enumerate(legend):
    if L[0] == None:
        ax_cb.text(0,-j,L[1],fontsize = fontsize, va='center')
    else:
        ax_cb.plot(0.2,-j+0.09,".",marker="s",ms=fontsize,color=L[0], markeredgewidth=.0)
        ax_cb.text(0.33, -j, L[1], fontsize=fontsize, va='center')
fig.savefig(data_folder + "/CLASSIFIER/RNA" + "/PLOT_expression-heatmap.png", dpi=300)
plt.show()