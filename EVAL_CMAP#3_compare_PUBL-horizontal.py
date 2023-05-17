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
from sklearn.decomposition import PCA
import cmapPy.clue_api_client.clue_api_client as cac
from cmapPy.clue_api_client.gene_queries import are_genes_in_api
import subprocess
import matplotlib as mpl
import matplotlib.gridspec as gridspec
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

#########################
data_folder = "*PAULA"
cluster = int(raw_input("... cluster k= "))
top = 5
drugs_selected = [0,1,2,3,4,5,13]
pathways_selected = [0,1,2,6,10]
investigate = "KU-0063794"
investigate_cluster = 0
remove_duplicates = True
combine_middle_clusters = False #bool(raw_input("... combine middle clusters? <bool> "))
#########################

if remove_duplicates == True:
    drugs = pd.read_excel(data_folder + "/CMAP/*RES_CMAP-filtered.xlsx")
    drugs = drugs.groupby(["name"], as_index=False).min()
    drugs.to_excel(data_folder + "/CMAP/*RES_CMAP-filtered_duplicates-removed.xlsx", index=False)
else:
    drugs = pd.read_excel(data_folder + "/CMAP/*RES_CMAP-filtered_duplicates-removed.xlsx")
pathways = pd.read_excel(data_folder + "/CMAP/*RES_CMAP-filtered-pathway.xlsx")
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")
samples = []
for sample in drugs.columns:
    if "U-" in sample:
        samples.append(sample[:-7])

### specific stages
sset = list(data_clinical["ID_set"].where(data_clinical["cluster k=%s"%cluster] == investigate_cluster).dropna())
sset = list(set(sset) & set(samples))
stages = []
for sample in sset:
    stages.append(list(data_clinical["type"].where(data_clinical["ID_set"] == sample).dropna())[0])
values = np.array(drugs[[x + "_paired" for x in sset]].where(drugs["name"] == investigate).dropna(thresh=1))[0]
print values
res = pd.DataFrame()
res["sample"] = sset
res["values"] = values
res["stages"] = stages
print res
res.to_excel(data_folder + "/CMAP/RES_%s_cluster-%s.xlsx"%(investigate, investigate_cluster))

clusters = []
if combine_middle_clusters == True:
    data_clinical["cluster k=%s" % cluster] = data_clinical["cluster k=%s" % cluster].replace([1,2,3,4],[1,1,1,2])
    max_cluster = 3
else:
    max_cluster = cluster

for c in range(max_cluster):
    sset = list(data_clinical["ID_set"].where(data_clinical["cluster k=%s"%cluster] == c).dropna())
    sset = list(set(sset) & set(samples))
    sset = [x + "_paired" for x in sset]
    print c, len(sset)
    #temp = drugs[["name"] + sset]
    for j,row in drugs.iterrows():
        drugs.at[j,"cluster %s"%c] = len(row[sset].dropna()) / float(len(row[sset]))
    #temp = temp.sort_values(by="filter", ascending=False)
    #temp.to_excel(data_folder + "/CMAP/TEMP_drugs_cluster-%s.xlsx"%c)
    clusters.append("cluster %s"%c)
for j,row in drugs.iterrows():
    samples_full = [x + "_paired" for x in samples]
    drugs.at[j,"all samples"] = len(row[samples_full].dropna()) / float(len(row[samples_full]))
drugs = drugs.sort_values(by="all samples",ascending=False)
drugs = drugs.reset_index(drop=True)
drugs.loc[:,["name"]+clusters+["all samples"]].to_excel(data_folder + "/CMAP/*RES_drugs.xlsx", index=False)
vs.style("science")
drugs = drugs.loc[:24,:]
#drugs = drugs.sort_values(by="all samples",ascending=True)
#drugs = drugs.reset_index(drop=True)

fig = plt.figure("PLOT",[12,3])
gs = gridspec.GridSpec(2, 5, width_ratios=[1,.02,.2,1,.02], height_ratios=[2,1], top=.4, bottom=.15, wspace=.1)
im = drugs.loc[:,clusters+["all samples"]]
print im
labels = drugs.loc[:,"name"]
for s in drugs_selected:
    labels[s] = labels[s] + " *"
ax = fig.add_subplot(gs[:,3])
#ax.set_title("Drugs", fontsize=10.)
ax.imshow(im.transpose(),cmap="RdBu_r", vmin=.0, aspect="auto")
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=8., rotation=45., ha="left")
ax.set_yticks(range(len(clusters)+1))
if combine_middle_clusters == True:
    clusters = ["PAULA I","PAULA II", "PAULA III" ]
else:
    clusters = ["PAULA I","PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III" ]
ax.set_yticklabels(clusters+["all samples"],rotation=0.,ha="right", va='center',fontsize=8.)
ax.grid(False)
ax.tick_params(axis="x", bottom=False, top=False, labelbottom=False, labeltop=True)
ax.tick_params(axis="y", right=False, left=False)
#for s in drugs_selected:
#    ax.get_xticklabels()[s].set_fontweight("bold")

cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=np.min(np.min(im))*100, vmax=np.max(np.max(im))*100)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax = fig.add_subplot(gs[:,4])
cb = fig.colorbar(sm,
             cax=ax, orientation='vertical', label='')
cb.set_ticks([0,40])
cb.set_label(label='Fraction of samples [%]',weight='regular', fontsize=7.)


clusters = []
for c in range(max_cluster):
    sset = list(data_clinical["ID_set"].where(data_clinical["cluster k=%s"%cluster] == c).dropna())
    sset = list(set(sset) & set(samples))
    sset = [x + "_paired" for x in sset]
    print c, len(sset)
    #temp = drugs[["name"] + sset]
    for j,row in pathways.iterrows():
        pathways.at[j,"cluster %s"%c] = len(row[sset].dropna()) / float(len(row[sset]))
    #temp = temp.sort_values(by="filter", ascending=False)
    #temp.to_excel(data_folder + "/CMAP/TEMP_drugs_cluster-%s.xlsx"%c)
    clusters.append("cluster %s"%c)
for j,row in pathways.iterrows():
    samples_full = [x + "_paired" for x in samples]
    pathways.at[j,"all samples"] = len(row[samples_full].dropna()) / float(len(row[samples_full]))
pathways = pathways.sort_values(by="all samples",ascending=False)
pathways = pathways.reset_index(drop=True)
pathways.loc[:,["name"]+clusters+["all samples"]].to_excel(data_folder + "/CMAP/*RES_pathways.xlsx", index=False)
pathways = pathways.loc[:24,:]
#pathways = pathways.sort_values(by="all samples",ascending=True)
#pathways = pathways.reset_index(drop=True)
vs.style("science")
im = pathways.loc[:,clusters+["all samples"]]
print im
labels = pathways.iloc[:,0]
for s in pathways_selected:
    labels[s] = labels[s]+ " *"
print labels

ax = fig.add_subplot(gs[:,0])
#ax.set_title("Pathways", fontsize=10.)
ax.imshow(im.transpose(),cmap="RdBu_r", vmin=.0, aspect="auto")
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,  fontsize=8., rotation=45., ha="left")
ax.set_yticks(range(len(clusters)+1))
if combine_middle_clusters == True:
    clusters = ["PAULA I","PAULA II", "PAULA III" ]
else:
    clusters = ["PAULA I","PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III" ]
ax.set_yticklabels(clusters+["all samples"],rotation=0.,ha="right", va='center',fontsize=8.)
ax.grid(False)
ax.tick_params(axis="x", bottom=False, top=False, labelbottom=False, labeltop=True)
ax.tick_params(axis="y", right=False, left=False)
#for s in pathways_selected:
#    ax.get_xticklabels()[s].set_fontweight("bold")

cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=np.min(np.min(im))*100, vmax=np.max(np.max(im))*100)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax = fig.add_subplot(gs[:,1])
cb = fig.colorbar(sm,
             cax=ax, orientation='vertical', label='Fraction of samples [%]')
cb.set_ticks([0,50])
cb.set_label(label='Fraction of samples [%]',weight='regular', fontsize=7.)

fig.savefig(data_folder + "/CMAP/PLOT_PUBL_heatmap_combine-II=%s_horizontal.png"%combine_middle_clusters,dpi=300, transparent=True)