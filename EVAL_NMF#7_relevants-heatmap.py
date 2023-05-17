import numpy as np
import pylab as plt
import sys, os, inspect
directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append('/Users/dressler/Documents/Core/_PythonModules')
import visualize as vs
import xml.etree.ElementTree as elt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines import statistics as lstat
import string
import seaborn as sns

vs.style('science')
colors = vs.colors_indexed

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)



### set basic parameters ###
data_folder = "*PAULA"
k = int(input("k = "))
mode = 2
percentile_chosen = 90
cluster_names = ["PAULA I", "PAULA II", "PAULA III", "PAULA IV", "PAULA V"]
############################



modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
mode = modes[mode]

basis = pd.read_excel(data_folder + "/CLUSTER/_files/PP_%s_k%s_basis.xlsx"%(mode,k))
scores = pd.read_excel(data_folder + "/CLUSTER/_files/PP_%s_k%s_scores.xlsx"%(mode,k)).iloc[:,0]
selects = pd.read_excel(data_folder + "/CLUSTER/_files/PP_%s_k%s_selects.xlsx"%(mode,k)).iloc[:,0]
protein_names = pd.read_excel(data_folder + "/CLUSTER/_files/*input_data.xlsx").loc[:,["accession"]]

g = sns.clustermap(basis, cmap = 'RdBu_r', cbar_pos=None, dendrogram_ratio=(0.4,.05), figsize = (6,14)) #coolwarm
g.gs.update(left=0.05, right=0.65, bottom=.07)
new_indices = g.dendrogram_row.reordered_ind
labely = []
labely_pos = []
percentiles_up = []
percentiles_low = []
relevants = []
for i, b in enumerate(basis.iloc[0,:]):
    percentiles_up.append(np.percentile(basis.iloc[:, i], 95.))
    percentiles_low.append(np.percentile(basis.iloc[:, i], 5.))
for y, l in enumerate(new_indices):
    name = protein_names.iloc[l, 0].split(" ")
    score = scores[l]
    select = selects[l]
    values = basis.iloc[l, :]
    #if (np.where(values > percentiles_up, True, False)).any() == True or (np.where(values < percentiles_low, True, False)).any() == True:
    if score > np.percentile(scores, percentile_chosen):
    #if select == True:
        relevants.append(protein_names.iloc[l, 0])
        if len(name) > 3:
            name = name[:4]
        name = string.join(name, " ")
        labely.append(name)
        labely_pos.append(y)
print labely
g.ax_heatmap.set_yticks(labely_pos)
g.ax_heatmap.set_yticklabels(labely, fontsize=3.5, rotation=45., ha='left', va='bottom')
g.ax_heatmap.set_xticklabels(cluster_names, fontsize=10., rotation=90., ha='right')
g.savefig(data_folder + "/CLUSTER/figures/PLOT_heatmap_k%s_%s-score-percentile.png"%(k, percentile_chosen),dpi=500)
relevants = pd.DataFrame(relevants)
relevants.to_excel(data_folder + "/CLUSTER/RES_cluster_relevants_%s-score-percentile.xlsx"%percentile_chosen, index=False)
plt.show()
