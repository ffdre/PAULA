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
import matplotlib as mpl
vs.style('science')
colors = vs.colors_indexed
markers = ["v","o"]

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)


### settings
data_folder = "*PAULA"
group_column = "binary_group"
set_column = "TMT set"
exclude = ["U-MIY-51_H"]
complete = False #True is not working
modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
colormap = "Spectral"
colors = ["grey", "black"]
zscore = True
mode_index = 2
subplots_x = 3
subplots_y = 3
show_names = False
show_legend = True


# allocate patients to groups
if complete == True: complete = "_COMPLETE-DATA"
else: complete = ""
data_file = data_folder + "/NORMICS_results/RES_data_%s%s.xlsx"%(modes[mode_index],complete)
protein_names = pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins%s.xlsx"%complete).iloc[:, 0] #_annotated
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical.xlsx")
try: os.mkdir(data_folder + "/PCA")
except: pass

### single plots

data = pd.read_excel(data_file)
if exclude != None:
    data = data.drop(columns = exclude)
#data = data.dropna(thresh=len(data.columns)*.9)
#data = data.reset_index()
data = 2.**data
patients = data.columns[:]
groups = {} #dict containing patients per group
group_list = []
sets = {} #dict containing set size
for group in np.sort(data_clinical.loc[:, group_column].unique()):
    groups[group] = []
    group_list.append(group)
print group_list

for patient in patients:
    for ID, group, Set in zip(data_clinical.loc[:, "ID_set"], data_clinical.loc[:, group_column], data_clinical.loc[:, set_column]):
        if ID == patient:
            groups[group].append(patient)
            sets[ID] = Set

for i, r in data.iterrows():
    mean = np.nanmean(data.iloc[i, :])
    #data.iloc[i, :] = data.iloc[i, :].fillna(mean) - np.nanmean(data.iloc[i, :])
    if zscore == True:
        data.iloc[i, :] = stats.zscore(data.iloc[i, :].fillna(mean))
    else:
        data.iloc[i, :] = data.iloc[i, :].fillna(mean)

# 2D
markersize = 3.
fig = plt.figure("PLOT",[12,12])
gs = mpl.gridspec.GridSpec(4,4,hspace=.4, wspace=.5, width_ratios = [1.5,2,1,1.5])
ax = fig.add_subplot(gs[0,0])
pca = PCA(n_components=2)
PC = pca.fit_transform(np.transpose(data))
varx, vary = pca.explained_variance_ratio_
#ax = fig.add_subplot(121)
ax.axhline(y=0,color="black",alpha=.3,linestyle="--")
ax.axvline(x=0,color="black",alpha=.3,linestyle="--")
print len(PC)
print sets
for i,sample in enumerate(PC):
    print sample
    for group in groups.keys():
        if data.columns[i] in groups[group]:
            color_index = group_list.index(group)
            color_index_set = sets[data.columns[i]]
            marker_index = group_list.index(group)
            if show_names == True:
                ax.annotate(sets[data.columns[i]], xy=[sample[0], sample[1]], xytext=(0, -1), ha='center', va='center',
                            textcoords='offset points', size = 7, weight='bold')

            ax.plot(sample[0], sample[1], ".", marker=markers[marker_index], markersize=markersize, markeredgewidth=.0,
                    alpha=.9, color=mpl.cm.get_cmap('Spectral')(color_index_set / float(max(data_clinical[set_column]))))
            ax.plot(sample[0], sample[1], ".", marker=markers[marker_index], markersize=markersize,
                    markeredgewidth=markersize / 5., alpha=1., \
                    markeredgecolor=colors[color_index], color="None")

ax.grid(False)

if show_legend == True:
    for j,group in enumerate(["Healthy", "Tumor"]):
        #color_index = group_list.index(group)
        ax.plot([], [], marker=markers[j], markersize=markersize, markeredgewidth=markersize/5., alpha=1., \
                    markeredgecolor=colors[j], color="None", label=group)

    ax.legend(fancybox=False)
ax.set_xlabel("Principal component 1 (%s %%)"%np.round(varx*100,2))
ax.set_ylabel("Principal component 2 (%s %%)"%np.round(vary*100,2))
fig.savefig(data_folder + "/PCA/PLOT_PCA-%s_NEW.png"%modes[mode_index], dpi=300, transparent=True)
plt.show()


sys.exit()
# 2D 2
fig = plt.figure("PCA %s 2 vs. 3" % data_folder, [8, 8])
fig.subplots_adjust(hspace=.4)
ax = fig.add_subplot(111)
pca = PCA(n_components=3)
PC = pca.fit_transform(np.transpose(data))
varx, vary, varz = pca.explained_variance_ratio_
#ax = fig.add_subplot(121)
ax.plot([0,0],[min(PC[:,2]),max(PC[:,2])],color="darkgrey")
ax.plot([min(PC[:,1]),max(PC[:,1])],[0,0],color="darkgrey")
print len(PC)
for i,sample in enumerate(PC):
    for group in groups.keys():
        if data.columns[i] in groups[group]:
            color_index = group_list.index(group)
            ax.plot(sample[1], sample[2], ".",markersize=10, alpha = .9, color = colors[color_index])
            if show_names == True:
                ax.annotate(sets[data.columns[i]], xy=[sample[1], sample[2]], xytext=(7, 0), ha='left', va='center',
                            textcoords='offset points')
ax.set_xlabel("Principal component 2 (%s %%)"%np.round(vary*100,2))
ax.set_ylabel("Principal component 3 (%s %%)"%np.round(varz*100,2))
ax.grid(False)
fig.savefig(data_folder + "/PCA/PLOT_PCA-%s_2+3.png"%modes[mode_index], dpi=300)

# 3D
fig = plt.figure("PCA %s" % data_folder, [8, 8])
fig.subplots_adjust(hspace=.4)
ax = fig.add_subplot(111, projection="3d")
pca = PCA(n_components=3)
PC = pca.fit_transform(np.transpose(data))
varx, vary, varz = pca.explained_variance_ratio_
#ax = fig.add_subplot(121)
ax.plot([0,0],[min(PC[:,1]),max(PC[:,1])],color="darkgrey")
ax.plot([min(PC[:,0]),max(PC[:,0])],[0,0],color="darkgrey")
ax.plot([0,0],[0,0],[min(PC[:,2]),max(PC[:,2])],color="darkgrey")
print len(PC)
for i,sample in enumerate(PC):
    for group in groups.keys():
        if data.columns[i] in groups[group]:
            color_index = group_list.index(group)
            ax.plot([sample[0]], [sample[1]], ".",zs=[sample[2]],markersize=10, alpha = .9, color = colors[color_index])
            if show_names == True:
                ax.text(sample[0], sample[1], sample[2], " %s"%sets[data.columns[i]], ha='left', va='center')
ax.set_xlabel("Principal component 1 (%s %%)"%np.round(varx*100,2))
ax.set_ylabel("Principal component 2 (%s %%)"%np.round(vary*100,2))
ax.set_zlabel("Principal component 3 (%s %%)"%np.round(varz*100,2))
ax.grid(False)
fig.savefig(data_folder + "/PCA/PLOT_PCA-%s_3D.png"%modes[mode_index], dpi=300)

plt.show()
sys.exit()
### comparison plot
for mode_index in range(7):
    data_file = data_folder + "/NORMICS_results/RES_data_%s.xlsx"%modes[mode_index]
    data = pd.read_excel(data_file)
    data = 2.**data
    patients = data.columns[:]
    groups = {} #dict containing patients per group
    group_list = []
    sets = {} #dict containing set size
    for group in data_clinical.iloc[:, 1].unique():
        groups[group] = []
        group_list.append(group)
    print group_list

    for patient in patients:
        for ID, group, Set in zip(data_clinical.iloc[:, 0], data_clinical.iloc[:, 1], data_clinical.iloc[:, 2]):
            if ID == patient:
                groups[group].append(patient)
                sets[ID] = Set
            elif data_folder == "TCGA" and ID in patient:
                groups[group].append(patient)

    for i, r in data.iterrows():
        mean = np.nanmean(data.iloc[i, :])
        #data.iloc[i, :] = data.iloc[i, :].fillna(mean) - np.nanmean(data.iloc[i, :])
        data.iloc[i, :] = stats.zscore(data.iloc[i, :].fillna(mean))

    fig = plt.figure("PCA %s compare normalization algorithms" % data_folder, [10, 10])
    fig.subplots_adjust(hspace=.4)
    fig.subplots_adjust(wspace=.33)
    ax = fig.add_subplot(subplots_y, subplots_x, mode_index+1)
    pca = PCA(n_components=2)
    PC = pca.fit_transform(np.transpose(data))
    varx, vary = pca.explained_variance_ratio_
    #ax = fig.add_subplot(121)
    ax.plot([0,0],[min(PC[:,1]),max(PC[:,1])],color="darkgrey")
    ax.plot([min(PC[:,0]),max(PC[:,0])],[0,0],color="darkgrey")
    print len(PC)
    for i,sample in enumerate(PC):
        for group in groups.keys():
            if data.columns[i] in groups[group]:
                color_index = group_list.index(group)
                ax.plot(sample[0], sample[1], ".",markersize=10, alpha = .9, color = colors[color_index])
                #ax.annotate(sets[data.columns[i]], xy=[sample[0], sample[1]], xytext=(7, 0), ha='left', va='center',
                 #           textcoords='offset points')
    ax.set_xlabel("Principal component 1 (%s %%)"%np.round(varx*100,2))
    ax.set_ylabel("Principal component 2 (%s %%)"%np.round(vary*100,2))
    ax.set_title(modes[mode_index])
    ax.grid(False)
fig.savefig(data_folder + "/PCA/PLOT_PCA-norm-compare.png", dpi=300)
plt.show()
sys.exit()


plt.show()
sys.exit()