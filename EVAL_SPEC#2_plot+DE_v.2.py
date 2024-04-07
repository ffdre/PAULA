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
import statistics_p3 as st
import string
import random
#from sklearn.decomposition import PCA
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import seaborn as sns
vs.style('science')
colors = vs.colors_indexed
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)


#######################
data_folder = "*PAULA"
parameter = "_UP_share"
parameter_mask = "_DN_share"
parameter_FC = "_UP_mean"
parameter_subplot1 = "_UP_share"
sort_by = "ALL"
name = "Upregulated targets"
top = 74
membrane = str(input("... membrane proteins only? <y/n> "))
in_manuscript = ["TMSB10","S100P", "S100A2", "TYMP", "CDK1", "TAGLN2",
                 "PLOD2", "MAL2", "LAMP2", "THEM6", "SLC9A1", "SLC14A1", "NECTIN4", "THBS1", "THBS2"]
#######################

def plot(ax,x,y,color,len_all):
    """Plot boxplot with swarmplot"""
    vp = ax.violinplot([y], positions=[x], showmeans=False, showmedians=False, showextrema=False, widths = .9)
    for param in vp['bodies']:
        param.set_alpha(0.25)
        param.set_facecolor(color)
    ax.boxplot([y], positions=[x], showfliers=False, patch_artist=True, widths=.45,
               boxprops=dict(facecolor=color, alpha=.25, linewidth=.5),  # , color=c),
               capprops=dict(color=color, alpha=.5),
               whiskerprops=dict(color=color, alpha=.5),
               flierprops=dict(color="black", markeredgecolor=None, markersize=5., marker="."),
               medianprops=dict(color="grey"))
    #ax.bar([x],[np.median(y)], color=color, width=.5, alpha=.25)
    if np.median(y) < np.mean(y):
        ax.plot([x-.2, x+.2],[np.mean(y), np.mean(y)], color=color, lw=2., alpha=.5)
    else:
        ax.plot([x-.2, x+.2],[np.mean(y), np.mean(y)], color=color, lw=2., alpha=.5)
        #ax.plot([x - .25, x + .25], [np.mean(y), np.mean(y)], color="white", lw=2., alpha=1.)
    sns.swarmplot(np.array([x] * len(y)), y, order=range(len_all), color=color, alpha=.5, ax=ax, size=3)
    ax.set_xticks(range(len_all))
    ax.set_xlim(-1,len_all)
    #ax.set_yscale('log')
    ax.grid(False)
    pass

sort_by += parameter
data = pd.read_excel(data_folder + "/SPECIFICITY/*RES_specificity_merged.xlsx")
if membrane == "y":
    data = data.where(data["membrane"] == 1).dropna(thresh=1)
groups_stage = ["ALL","pTa low", "pTa high", "pTis", "pT1a/b/G1/G2", "pT1c/G3", "pT2", ">pT2"] #,"pN1"
groups_cluster = ["PAULA I", "PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III"]
groups = groups_stage + groups_cluster

fig = plt.figure("PLOT", [12,9])
fig.subplots_adjust(hspace=.0, wspace=.2, bottom=.15)

#data["sort"] = data[["%s%s"%(x, parameter) for x in groups[1:]]].max(axis=1, skipna=True)
data["sort"] = data["%s%s"%("ALL", parameter)]
data = data.sort_values(by="sort", ascending=False)
data = data.reset_index(drop=True)



gs = gridspec.GridSpec(32, 7, height_ratios=list(np.array([1.]*30)) + [.4,.4], width_ratios=[1,7,5,.25,3,.5,3.5])
## data
im = data.loc[:top,["%s%s"%(x, parameter) for x in groups]]#.replace(np.nan,0.)
im_mask = data.loc[:top,["%s%s"%(x, parameter_mask) for x in groups]].values

data.loc[:,["Gene Symbol"]+["%s%s"%(x, parameter) for x in groups]].\
    to_excel(data_folder + "/SPECIFICITY/PLOT_values_membrane=%s.xlsx"%membrane, index=False)
data.loc[:,["Gene Symbol"]+["%s%s"%(x, parameter_mask) for x in groups]].\
    to_excel(data_folder + "/SPECIFICITY/PLOT_values-mask_membrane=%s.xlsx"%membrane, index=False)

im_fused = []
for k, ii in enumerate(im.values):
    temp = []
    for l, iii in enumerate(im.values[k]):
        temp.append(im.values[k,l])
        temp.append(-im_mask[k,l])
    im_fused.append(temp)
im_fused = np.array(im_fused)


## Subplot variance
print(len(data.iloc[:,0]))
pTNM = data.loc[:top,["%s%s"%(x, parameter_subplot1) for x in groups_stage[1:]]].var(axis=1).values * 7./7.
paula = data.loc[:top,["%s%s"%(x, parameter_subplot1) for x in groups_cluster[:]]].var(axis=1).values * 7./7.
st.deltaquant(pTNM,paula,names=["pTNM","PAULA"])

ax = fig.add_subplot(gs[:8,-3:-2])
#print im.loc[:,["%s%s"%(x, parameter) for x in groups_stage[1:]]].var(axis=1)
#for a,b in zip(pTNM,paula):
#    ax.plot([0,1],[a,b],alpha=.1, color="black", lw=.5)
plot(ax,0,pTNM,colors[0],2)
plot(ax,1,paula,colors[1],2)
ax.set_ylabel("Variance [1]")
ax.set_xticks([0,1])
ax.set_xticklabels(["pTNM","PAULA"],rotation=45.,ha="right")
ax.set_xlim(-.6,1.6)
ax.text(.05, .95, "p = 0.006", transform=ax.transAxes, ha='left', va='top')
ax.grid(False)

## Subplot FC
print(len(data.iloc[:,0]))
pTNM = data.loc[:top,["%s%s"%(x, parameter_FC) for x in groups_stage[1:]]].mean(axis=1, skipna=True).values * 7./7.
paula = data.loc[:top,["%s%s"%(x, parameter_FC) for x in groups_cluster[:]]].mean(axis=1, skipna=True).values * 7./7.
st.deltaquant(pTNM,paula,names=["pTNM_FC","PAULA_FC"])

ax = fig.add_subplot(gs[12:20,-3:-2])
#print im.loc[:,["%s%s"%(x, parameter) for x in groups_stage[1:]]].var(axis=1)
plot(ax,0,pTNM,colors[0],2)
plot(ax,1,paula,colors[1],2)
ax.set_ylabel("Mean FC [log$_2$]")
ax.set_xticks([0,1])
ax.set_xlim(-.6,1.6)
ax.set_xticklabels(["pTNM","PAULA"],rotation=45.,ha="right")
ax.text(.05, .95, "p = n.s.", transform=ax.transAxes, ha='left', va='top')
ax.grid(False)

#ax = fig.add_subplot(gs[5:10,-1])
#print im.loc[:,["%s%s"%(x, parameter) for x in groups_stage[1:]]].var(axis=1)
#ax.boxplot([im.loc[:,["%s%s"%(x, parameter) for x in groups_stage[1:]]].values.flatten(),
#           im.loc[:,["%s%s"%(x, parameter) for x in groups_cluster[:]]].values.flatten()])
#ax.grid(False)

#print("... variance pTNM: %s"%np.mean(im.loc[:,["%s%s"%(x, parameter) for x in groups_stage[1:]]].var()))
#print("... variance PAULA: %s"%np.var(.values.flatten()))

labels = data.loc[:top,"Gene Symbol"]
labels = [x+" *" if x in in_manuscript else x for x in labels ]
vmax = np.max(im.values)
reds = plt.cm.get_cmap('Reds', 512)
#Reds_new = mpl.colors.ListedColormap(reds(np.linspace(0., .6, 256)))
cdict = {'red':   ((.0,  1., 1.), #(0, 114, 197)
                   (1.,  .0, .0)),
         'green':  ((.0,  .0, .0),
                    (.0, .0, .0)),
         'blue':  ((.0,  .0, .0), #(220, 68, 56)
                   (1.,  1., 1.))}
Reds_new = mpl.colors.LinearSegmentedColormap.from_list('Reds_new',[colors[0],"#e3bfbc", colors[1]])
Traffic_light = mpl.colors.LinearSegmentedColormap.from_list('Traffic_light',[vs.colors_cluster[4],"black", vs.colors_cluster[1]])
Traffic_light = mpl.colors.LinearSegmentedColormap.from_list('Traffic_light',["#fc3103","black", "#03fc62"])

cmap = Traffic_light
## Subplot ALL
ax = fig.add_subplot(gs[:-2,0])
axC = fig.add_subplot(gs[-1,0])
axC.tick_params(axis="x", bottom=True, top=False, labelbottom=False, labeltop=False)
axC.tick_params(axis="y", right=False, left=False, labelleft=False)
axC.axis('off')
ax.set_title("All")
ax.imshow(im_fused[:,0:2], cmap=cmap, vmin=-.5,vmax=+.5, aspect="auto")
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels, fontsize=7.)
ax.set_xticks(np.arange(len(groups[0:1]))*2+1)
ax.set_xticklabels(["All samples"],rotation=45.,ha="right")
ax.grid(False)
ax.tick_params(axis="x", bottom=False, top=False, labelbottom=True, labeltop=False)
ax.tick_params(axis="y", right=False, left=False)

## Subplot pTNM
ax = fig.add_subplot(gs[:-2,1])
ax.set_title("pTNM")

axC = fig.add_subplot(gs[-1,1])
axC.imshow([np.arange(len(groups[1:8]))],cmap="Spectral_r", aspect="auto", vmax=len(groups[1:8]))
axC.tick_params(axis="x", bottom=False, top=False, labelbottom=True, labeltop=False)
axC.tick_params(axis="y", right=False, left=False, labelleft=False)
axC.grid(False)

ax.imshow(im_fused[:,2:16], cmap=cmap, vmin=-.5,vmax=+.5, aspect="auto")
axC.set_xticks(np.arange(len(groups[1:8])))
axC.set_xticklabels(groups[1:8],rotation=45.,ha="right")
ax.grid(False)
ax.tick_params(axis="x", bottom=False, top=False, labelbottom=False, labeltop=False)
ax.tick_params(axis="y", right=False, left=False, labelleft=False)

## Subplot PAULA
ax = fig.add_subplot(gs[:-2,2])
ax.set_title("PAULA")

axC = fig.add_subplot(gs[-1,2])
axC.imshow([np.arange(len(groups[8:]))],cmap=mpl.colors.ListedColormap(vs.colors_cluster), aspect="auto")
axC.tick_params(axis="x", bottom=False, top=False, labelbottom=True, labeltop=False)
axC.tick_params(axis="y", right=False, left=False, labelleft=False)
axC.grid(False)

ax.imshow(im_fused[:,16:], cmap=cmap, vmin=-.5,vmax=+.5, aspect="auto")
axC.set_xticks(np.arange(len(groups[8:])))
axC.set_xticklabels(groups[8:],rotation=45.,ha="right")
ax.grid(False)
ax.tick_params(axis="x", bottom=False, top=False, labelbottom=False, labeltop=False)
ax.tick_params(axis="y", right=False, left=False, labelleft=False)

#cmap = Traffic_light #mpl.cm.Reds #RdBu_r
norm = mpl.colors.Normalize(vmin=-50., vmax=+50.)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax = fig.add_subplot(gs[-7:-1,3:6])
ax.grid(False)
ax.axis('off')
cbar = fig.colorbar(sm, fraction=1.,
             ax=ax, orientation='vertical')
cbar.set_ticks([-50,0,50])
cbar.ax.set_yticklabels([50,0,50])
cbar.ax.tick_params(size=0) #Remove ticks
cbar.outline.set_visible(False) #Remove outline
cbar.set_label('Fraction of samples [%] with\nunderexpression | overexpression', fontsize=8)

fig.savefig(data_folder + "/SPECIFICITY/PLOT_PUBL_DOUBLE_%s_membrane=%s_v.2.png"%(parameter,membrane), dpi=300, transparent=True)
plt.show()
