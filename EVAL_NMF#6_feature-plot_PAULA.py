
from collections import defaultdict, Counter
import urllib
import sys, os, inspect
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch
directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
sys.path.append(file_directory + '/_modules')
import nimfa
import visualize as vs
import pandas as pd
from scipy import stats
import scipy
import matplotlib as mpl
vs.style('science')
colors = vs.colors_indexed

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)



### set basic parameters ###
data_folder = "*PAULA"
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")
mode_index = 2
cluster = 5

fontsize = 10
#############################

legend = []
features = []
cmaps = ['tab10','Spectral'] #has to be of length cohort_features Set1, tab10, tab20  # "Set1", "Spectral" for Lung
modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedNORMICSmed-log2_PM",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
mode = modes[mode_index]
samples = np.load(data_folder + "/CLUSTER/_files" + "/PP_patients.npy")
data_clinical = data_clinical.sort_values(by = "image_order k=%s"%cluster)
data_clinical = data_clinical.where(data_clinical["ID_set"].isin(samples)).dropna(thresh=1)

vs.style('science')
colors = vs.colors_indexed

fig = plt.figure("features",[13.9, 12])
fig.tight_layout()
fig.subplots_adjust(left=.1, right=.9)
GS = gridspec.GridSpec(38, 2, figure=fig, wspace=.0, hspace=0.2, #0.018
                              width_ratios=[0.4, 1])

ax_cb = fig.add_subplot(GS[:,0])
ax_cb.axis('off')

def add_feature(AX, data, name, cmap, names, silent=False):
    global legend, features
    elements = np.sort(data.unique())
    if silent == False:
        legend.append((None, name))
        for n, e in zip(names, elements):
            legend.append((mpl.cm.get_cmap(cmap)((e - np.min(elements)) / float(np.max(elements) - np.min(elements))), n))
        legend.append((None, ""))
    AX.imshow([data], cmap=cmap, aspect="auto")
    AX.set_yticks([0])
    AX.set_yticklabels([name], fontsize=fontsize)
    features.append(AX)
    pass

def add_feature_binary(AX, data, name, colors, names, silent=False):
    global legend, features
    elements = np.sort(data.unique())
    if silent == False:
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
data = data_clinical["TMT_set_number"]
name = "TMT set"
names = range(len(data.unique()))
cmap = "Spectral"  # "Spectral" Set1
add_feature(ax, data, name, cmap, names, silent=True)

ax = fig.add_subplot(GS[1, 1])
data = data_clinical["TMT_set_low"]
name = "TMT low"
cmap = [vs.colors["blue"], vs.colors["darkestblue"]]
names = ["no", "yes"]
add_feature_binary(ax, data, name, cmap, names)

ax = fig.add_subplot(GS[2, 1])
data = data_clinical["Stage"]
name = "pTNM stage"
cmap = "Spectral_r"
names = ["pTa low", "pTa high", "pTis", "pT1a/b/G1/G2", "pT1c/G3", "pT2", ">pT2","pN1"]
add_feature(ax, data, name, cmap, names)

ax = fig.add_subplot(GS[3, 1])
data = data_clinical["cluster k=5"]
name = "Proteomic cluster"
cmap = vs.colors_cluster
names = ["PAULA I", "PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III"]
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

# colorbar heatmap
cmap = mpl.cm.RdBu
norm = mpl.colors.Normalize(vmin=0., vmax=1.)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax = fig.add_subplot(4,120,363)
cb = fig.colorbar(sm,cax=ax, orientation='vertical') #label='fold change [log2]'
cb.ax.tick_params(labelsize=8)
cb.outline.set_linewidth(.0)
ax.set_ylabel("inter-run distance [1]", fontsize=fontsize)

fig.savefig(data_folder + "/CLUSTER/figures" + "/PLOT_features.png", dpi=300)
plt.show()