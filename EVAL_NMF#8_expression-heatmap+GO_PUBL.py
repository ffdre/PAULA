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
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy.cluster.hierarchy import dendrogram, linkage
vs.style('science')
colors = vs.colors_indexed

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)


### settings ############################
data_folder = "*PAULA"
reverse_log2 = False
zscore = True
fill_na = False
clusters = int(input("k = "))
percentile = int(input("percentile = "))
mode_index = 2
num_features = 5
len_labels = 40
max_heatmap = 2

fontsize = 10
fontsize_heatmap = 5
legend = []
features = []
modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
mode = modes[mode_index]
#########################################


def add_feature(AX, data, name, cmap, names, silent=False, vertical=False, plot=True):
    global legend, features
    elements = np.sort(data.unique())
    if silent == False:
        legend.append((None, name))
        for n, e in zip(names, elements):
            legend.append((mpl.cm.get_cmap(cmap)((e - np.min(elements)) / float(np.max(elements) - np.min(elements))), n))
        legend.append((None, ""))
    if plot == True:
        if vertical == False:
            AX.imshow([data], cmap=cmap, aspect="auto")
            AX.set_yticks([0])
            AX.set_yticklabels([name], fontsize=fontsize)
            AX.set_xticklabels([])
        else:
            AX.imshow(np.transpose(np.array([data])), cmap=cmap, aspect="auto")
            AX.set_xticks([0])
            AX.set_xticklabels([name], fontsize=fontsize, rotation=45., ha='right')
            AX.set_yticklabels([])
        features.append(AX)
    pass

def add_feature_binary(AX, data, name, colors, names, silent=False, vertical=False, plot=True):
    global legend, features
    elements = np.sort(data.unique())
    if silent == False:
        legend.append((None, name))
        for n, e in zip(names, elements):
            legend.append((colors[int(e)], n))
        legend.append((None, ""))
    if plot == True:
        if vertical == False:
            AX.imshow([data], cmap=mpl.colors.ListedColormap(colors), aspect="auto")
            AX.set_yticks([0])
            AX.set_yticklabels([name], fontsize=fontsize)
            AX.set_xticklabels([])
        else:
            AX.imshow(np.transpose(np.array([data])), cmap=mpl.colors.ListedColormap(colors), aspect="auto")
            AX.set_xticks([0])
            AX.set_xticklabels([name], fontsize=fontsize, rotation=45., ha='right')
            AX.set_yticklabels([])
        features.append(AX)
    pass

data = pd.read_excel(data_folder + "/CLUSTER/_files/*input_data.xlsx")
data = pd.merge(data, pd.read_excel(data_folder + "/*RES_data_functional+type.xlsx")[["Accession","Gene Symbol","Biological Process","membrane"]],
                how='left', left_on="accession", right_on="Accession")
sample_columns = ~data.columns.isin(['sort', 'accession','Accession','Gene Symbol','Biological Process','membrane'])
## reverse log2
if reverse_log2 == True:
    data.loc[:,sample_columns] = 2. ** data.loc[:,sample_columns]

## fill NA, prepare sorting, transform to zscore
for i, row in data.iterrows():
    mean = np.nanmean(data.iloc[i, sample_columns])
    SD = np.nanstd(data.iloc[i, sample_columns], ddof=0)
    if zscore == True:
        data.iloc[i, sample_columns] = (data.iloc[i, sample_columns] - mean)/float(SD)
        if fill_na == True:
            data.iloc[i, sample_columns] = data.iloc[i, sample_columns].fillna(0)
    else:
        if fill_na == True:
            data.iloc[i, sample_columns] = data.iloc[i, sample_columns].fillna(mean)

basis = pd.read_excel(data_folder + "/CLUSTER/_files/PP_%s_k%s_basis.xlsx" % (mode, clusters))
#print basis.columns
sorters = []
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")
data_clinical = data_clinical.where(data_clinical["ID_set"].isin(data.columns)).dropna(thresh=1)
data_clinical = data_clinical.sort_values(by=["cluster k=%s"%clusters], ascending=True)
sorted_samples = data_clinical["ID_set"]
cluster_samples = data_clinical["cluster k=%s"%clusters].values
print(data_clinical["cluster k=%s"%clusters].unique())
for k in data_clinical["cluster k=%s"%clusters].unique():
    if math.isnan(k) == False:
        data["cluster %s"%k] = basis.iloc[:,int(k)] * data[data_clinical["ID_set"].where(\
            data_clinical["cluster k=%s"%clusters] == k).dropna()].median(axis=1)
        sorters.append("cluster %s"%k)

relevants = pd.read_excel(data_folder + "/CLUSTER/RES_cluster_relevants_%s-score-percentile.xlsx"%percentile).iloc[:,0]
data = data.where(data["accession"].isin(relevants)).dropna(thresh=1)

### get frequencies of GO terms
GOs = []
for GO in data["Biological Process"].values:
    terms = str(GO).split(";")
    for t in terms: GOs.append(t)
GOs = np.array(GOs)

GOs = np.unique(GOs, return_counts=True)
print(GOs)
GOs = pd.DataFrame(GOs).transpose()
GOs = GOs.sort_values(by=GOs.columns[1],ascending=True)
GO_terms = GOs.iloc[:,0].values
GOs.to_excel(data_folder + "/CLUSTER/figures" + "/RES_expression-heatmap_%spercentile-relevants_GO-list.xlsx"%percentile,index=False)
### end

#data = data.sort_values(by=sorters, ascending=False)
data = data.reset_index(drop=True)
data_vis = data.loc[:, sample_columns]
data_vis = data[sorted_samples]
#data_vis = data_vis.set_index(data["Gene Symbol"].values)

def clean_axis(ax, y=False):
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().grid(False)
    if y == False: ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values(): sp.set_visible(False)


fig = plt.figure("expression heatmap",[15, 13])
fig.tight_layout()
fig.subplots_adjust(left=.1, right=.9)

Y = linkage(data_vis.groupby(by=cluster_samples, axis=1).mean(), method='average') #.groupby(by=cluster_samples, axis=1).mean()
width_for_GOs = .15
number_GOs = 12
GS2 = gridspec.GridSpec(38, 2 + number_GOs, figure=fig, wspace=.05, hspace=0.2, #0.018
                              width_ratios=[0.4-width_for_GOs]+list(np.array([width_for_GOs/float(number_GOs)]*number_GOs)) + [1.2])
ax = fig.add_subplot(GS2[num_features:,1])
ax.grid(False)
clean_axis(ax, y=True)
#ax.axis("off")
#ax.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True)
denD = dendrogram(Y, orientation='left', link_color_func=lambda k: 'black', labels=data["Gene Symbol"].values,no_labels=False, no_plot=True)
reorder = denD['leaves']
data_vis = data_vis.iloc[reorder,:]

GO_terms = [GO for GO in GO_terms if GO != 'nan' and 'other bio' not in GO]


for f,GO in enumerate(GO_terms):
    ax = fig.add_subplot(GS2[num_features:, f+1])
    d = []
    for da in data.iloc[reorder,:]["Biological Process"].values:
        if GO in str(da): d.append(1.)
        else: d.append(0.)
    name = ""
    cmap = ["lightgrey", mpl.cm.get_cmap("Spectral")(f/float(len(GO_terms)-1))]
    names = ["no", "yes"]
    add_feature_binary(ax, pd.Series(d), name, cmap, names, silent=True, vertical=True)

GS = gridspec.GridSpec(38, 2, figure=fig, wspace=.0, hspace=0.2, #0.018
                              width_ratios=[0.4, 1])
ax = fig.add_subplot(GS[num_features:,1])
if max_heatmap == False:
    vmin = -np.max([abs(np.min(data_vis)), abs(np.max(data_vis))])
    vmax = np.max([abs(np.min(data_vis)), abs(np.max(data_vis))])
else:
    vmin, vmax = -max_heatmap, max_heatmap
cmap = mpl.cm.RdBu_r
cmap.set_bad('grey',1.)
ax.imshow(data_vis, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
ax.grid(False)
ax.set_yticks(range(len(data["Gene Symbol"])))
labels = []
for label in data["Gene Symbol"][reorder]:
    label = str(label)
    name = label.split(";")[0]
    if len(name) > len_labels:
        name = name[:len_labels] + "..."
    labels.append(name)
ax.set_yticklabels(labels, fontsize=fontsize_heatmap)
#ax.set_yticklabels(labels)
ax.tick_params(left=False, right=False, top=False, bottom=False, labelbottom=False, labelleft=False, labelright=True)
ax.grid(False)
ax.set_frame_on(False)
#features.append(ax)

fig2 = plt.figure(2)
ax = fig2.add_subplot(GS[num_features:,1])
#vmin = -np.max([abs(np.min(data_vis)), abs(np.max(data_vis))])
#vmax = np.max([abs(np.min(data_vis)), abs(np.max(data_vis))])
ax.imshow(data_vis.groupby(by=cluster_samples, axis=1).mean(), cmap=cmap, aspect="auto") #.groupby(by=cluster_samples, axis=1).mean()
data_vis_samples = data_vis.copy()
data_vis_samples["proteins"]=labels
data_vis_samples.to_excel(data_folder + "/CLUSTER/figures/RES_data_visualized_samples.xlsx",index=False)
data_vis = data_vis.groupby(by=cluster_samples, axis=1).mean()
data_vis["proteins"]=labels
data_vis.to_excel(data_folder + "/CLUSTER/figures/RES_data_visualized.xlsx",index=False)
ax.grid(False)

# colorbar heatmap
cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax = fig.add_subplot(4,120,363)
cb = fig.colorbar(sm,cax=ax, orientation='vertical') #label='fold change [log2]'
cb.ax.tick_params(labelsize=8)
cb.outline.set_linewidth(.0)
ax.set_ylabel("log$_2$ z-score [1]", fontsize=fontsize)

## Features
ax_cb = fig.add_subplot(GS[:,0])
ax_cb.axis('off')


### FEATURES ###
ax = fig.add_subplot(GS[0, 1])
data = data_clinical["TMT_set_number"]
name = "TMT set"
names = range(len(data.unique()))
cmap = "Spectral"  # "Spectral" Set1
add_feature(ax, data, name, cmap, names, silent=True)

ax = fig.add_subplot(GS[1, 1])
data = data_clinical["TMT_set_low"]
name = "TMT set type"
cmap = [vs.colors["darkestblue"], vs.colors["blue"]]
names = ["high", "low"]
add_feature_binary(ax, data, name, cmap, names)

ax = fig.add_subplot(GS[2, 1])
data = data_clinical["Stage"]
name = "pTNM stage"
cmap = "Spectral_r"
names = ["pTa low", "pTa high", "pTis", "pT1a/b/G1/G2", "pT1c/G3", "pT2", ">pT2","pN1"]
add_feature(ax, data, name, cmap, names)

ax = fig.add_subplot(GS[3, 1])
data = data_clinical["patient_ID"]
print("... # of patients: %s"%len(data_clinical["patient_ID"].unique()))
name = "Patient"
cmap = "Spectral_r"
names = range(len(data.unique()))
data = data.replace(data.unique(), names)
add_feature(ax, data, name, cmap, names, silent=True)

ax = fig.add_subplot(GS[4, 1])
data = data_clinical["cluster k=%s"%clusters]
name = "Proteomic cluster"
cmap = vs.colors_cluster
names = ["PAULA I", "PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III"]
add_feature_binary(ax, data, name, cmap, names)

name = "GO terms (ordinate)"
names = []
cs = 20
for g in GO_terms:
    g = g.replace(" OR ","/")
    g = g.replace(" and ","/")
    g = g.replace("organization","organiz.")
    g = g.replace("metabolism","metab.")
    g = g.replace("metabolic","metab.")
    if len(g) > cs and g[cs-1] != " ":
        names.append(g[:cs]+"...")
    elif len(g) > cs:
        names.append(g[:cs-2]+"...")
    else: names.append(g)
cmap = "Spectral"  # "Spectral" Set1
add_feature(ax, pd.Series(range(len(names))), name, cmap, names, silent=False, plot=False)

### rest ###
for ax in features:
    ax.grid(False)
    ax.set_frame_on(False)
    ax.tick_params(left=False, right=False, top=False, bottom=False)



ax_cb.set_ylim(-50,.4)
ax_cb.set_xlim(0,3)
for j,L in enumerate(legend):
    if L[0] == None:
        ax_cb.text(0,-j,L[1],fontsize = fontsize, va='center')
    else:
        ax_cb.plot(0.2,-j+0.09,".",marker="s",ms=fontsize,color=L[0], markeredgewidth=.0)
        ax_cb.text(0.33, -j, L[1], fontsize=fontsize, va='center')
fig.savefig(data_folder + "/CLUSTER/figures" + "/PLOT_expression-heatmap_%spercentile-relevants_GO.png"%percentile, dpi=300, transparent=True)
plt.show()
