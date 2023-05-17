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
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

#########################
data_folder = "*PAULA"
extract_first = True
#########################

data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical.xlsx")

if extract_first == True:
    counter = 0
    for file in sorted(os.listdir(os.getcwd()+"/%s/CMAP/data"%data_folder)):
        if file.endswith(".gct"):

            run = file.split("_")[1]
            #print file[2:]
            print(">>> loading batch #%s"%run)
            data_temp = pd.read_csv("%s/CMAP/data/"%data_folder+file, header=[0,1,2], delimiter="\t", skiprows=[0,1,5,6,7,8])#, index_col=concatenation_column) header=[4,3,2],
            data_temp = data_temp.where(data_temp[("type","na","na")] == "CP").dropna(subset=[("type","na","na")])
            columns_to_keep, columns_new_names = [], []
            for column in data_temp.columns:
                if column[1] == "summary" or column[0] in ["id", "name", "description", "target"]:
                    columns_to_keep.append(column)
                    if column[1] != "summary":
                        columns_new_names.append(column[0])
                    else:
                        columns_new_names.append(column[2])
            data_temp = data_temp[columns_to_keep]
            data_temp.columns = columns_new_names
            if counter == 0:
                data = data_temp
            else:
                data = pd.merge(data,data_temp, how="outer", on=["id", "name", "description", "target"])
            counter += 1
            #print data_temp
            #data.to_excel("%s/CMAP/data/"%data_folder+"check.xlsx")
            #sys.exit()
    #print data
    data.to_excel("%s/CMAP/"%data_folder+"*RES_CMAP.xlsx", index=False)
else:
    data = pd.read_excel("%s/CMAP/"%data_folder+"*RES_CMAP.xlsx")

samples = []
samples_ID = []
for column in data.columns:
    if "_paired" in column:
        samples.append(column)
        samples_ID.append(column[:-7])
data_masked = data[samples].where(data[samples]<=-90.)
data_masked[["id", "name", "description", "target"]] = data[["id", "name", "description", "target"]]
data = data_masked.dropna(thresh=1,subset=samples)
print data
data.to_excel("%s/CMAP/"%data_folder+"*RES_CMAP-filtered.xlsx", index=False)

pathways_not_clean = list(data["description"].unique())
pathways = []
for pathway_not_clean in pathways_not_clean:
    for pathway in pathway_not_clean.split(", "):
        if pathway not in pathways:
            pathways.append(pathway)

data_pathway = pd.DataFrame(index=pathways, columns = samples)

for pathway in data_pathway.index:
    temp = []
    indices = list(data.index[data["description"] == pathway])
    #print indices
    for sample in samples:
        count = len(data.loc[indices,sample].dropna())
        if count > 0:
            temp.append(count)
        else:
            temp.append(np.nan)
    #print temp
    #temp = pd.Series(temp, columns = samples)
    data_pathway.loc[pathway] = temp

print data_pathway
data_pathway.to_excel("%s/CMAP/"%data_folder+"*RES_CMAP-filtered-pathway.xlsx", index=True)




groups = ["healthy","pTa low", "pTa high", "pTis", "pT1a/b G2", "pT1c G3", "pT2", ">pT2", "pN1"]
# Stage ID_set


group_indices = []
for x, stage in enumerate(data_clinical["Stage"].sort_values(ascending=True).unique()):
    relevant_samples = list(data_clinical.where(data_clinical["Stage"] == stage).loc[:,"ID_set"].dropna())
    actual_samples = []
    for sample in relevant_samples:
        if sample in samples_ID:
            actual_samples.append(sample+"_paired")
    temp = []
    for j, row in data_pathway.iterrows():
        count = len(list(row[actual_samples].dropna()))
        temp.append(count/float(len(actual_samples)))
    group_indices.append("%s n=%s"%(groups[stage], len(actual_samples)))
    data_pathway["%s n=%s"%(groups[stage], len(actual_samples))] = temp
data_pathway.to_excel("%s/CMAP/" % data_folder + "*RES_CMAP-filtered-pathway.xlsx", index=True)

vs.style('science', n=16, colormap = "Spectral")
colors = vs.colors_indexed
fig = plt.figure("Pathway profiles across stages",[7,7])
for i, group_index in enumerate(group_indices):

    ax = fig.add_subplot(3,3,i+1)
    ax.set_title(group_index, fontsize=10, fontweight="regular")
    counter =0
    labels = []
    for index,row in data_pathway.iterrows():
        if max(row[group_indices]) > 0.3:
            ax.bar(counter, .6, color="lightgrey", linewidth=0.)
            ax.bar(counter, row[group_index], color = colors[counter], linewidth=0.)
            #ax.plot(group_indices, row[group_indices], color = colors[counter], alpha=.75)
            labels.append(index)
            counter += 1
    ax.grid(False)
    ax.set_ylim(0,0.6)
    ax.tick_params(top='off', bottom='off', left='on', right='off', labelleft='on', labelbottom='off')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    #else:
    #    ax.plot(group_indices, row[group_indices], label="_nolegend_",color = "grey", alpha=.25)
#ax.legend(fancybox=False)
#ax.xticklabels(rotation=45, ha="right")


fig = plt.figure("Legend",[.5,4])
for i,item in enumerate(labels):
    ax = fig.add_subplot(len(labels), 2, i * 2 + 2)
    ax.imshow([[i/float(len(labels))]], cmap="Spectral", vmin=0., vmax=1.)
    print i / float(len(labels))
    ax.yaxis.set_ticks(range(1))
    ax.yaxis.set_ticklabels([item])
    ax.xaxis.set_ticklabels([])
    ax.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')
    ax.grid(False)
plt.show()