
from collections import defaultdict, Counter
import urllib
import sys, os, inspect
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch
import pandas as pd
from scipy import stats
import scipy
import nimfa
from scipy.spatial.distance import squareform
sys.path.append('/Users/dressler/Documents/Core/_PythonModules')
import visualize as vs
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
vs.style('science')
colors = vs.colors_indexed
plt.close('all')

### set basic parameters ###
cluster_quantile_upper = .0
cluster_quantile_lower = .25
data_folder = "*PAULA"
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical.xlsx", index_col = "ID_set")
by_longlist = True
tumor_only = False
############################

### set adv. parameters ####
normed_to_healthy = False
complete_data = "" #"_COMPLETE-DATA"
mode_index = 2
zscore = False
log2 = True
cohort_features = ["Stage", "TMT_set_number"]
height_ratios = [0.025,0.025,0.9]
cmaps = ['Spectral','Spectral'] #has to be of length cohort_features Set1, tab10, tab20  # "Set1", "Spectral" for Lung
cmap_low = [None, None]
cmap_up = [None, None]
############################

## prepare
modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedNORMICSmed-log2_PM",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
mode = modes[mode_index]
try: os.mkdir(data_folder + "/CLUSTER")
except: pass
try: os.mkdir(data_folder + "/CLUSTER/figures")
except: pass
try: os.mkdir(data_folder + "/CLUSTER/_files")
except: pass

## write protocol
with open(data_folder + "/CLUSTER/_files/_protocol.txt","w") as f:
    for parameter, value in zip([cluster_quantile_upper, cluster_quantile_lower, data_folder, by_longlist, \
                                 tumor_only, normed_to_healthy, complete_data, modes[mode_index], zscore, \
                                 log2, cohort_features, height_ratios, cmaps, cmap_low, cmap_up], \
                                ["cluster_quantile_upper", "cluster_quantile_lower", "data_folder", "by_longlist", \
                                 "tumor_only", "normed_to_healthy", "complete_data", "mode", "zscore", \
                                 "log2", "cohort_features", "height_ratios", "cmaps", "cmap_low", "cmap_up"]):
        f.write("%s\t%s\n"%(parameter,value))

## load data
if normed_to_healthy == True:
    data = pd.read_excel(data_folder + "/NORMICS_results/RES_data_%s_PAIRED_%s.xlsx" %(modes[mode_index],complete_data))
else:
    data = pd.read_excel(data_folder + "/NORMICS_results/RES_data_%s%s.xlsx"%(modes[mode_index], complete_data))
data_longlist = pd.read_excel(data_folder + "/NORMICS_results/RES_longlist.xlsx")
data_protein_names = pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins%s.xlsx"%complete_data)
num_proteins = len(data_protein_names.iloc[:, 0])

## filter data
if tumor_only == True:
    actual_columns = []
    for column in data.columns:
        if "_H" not in column:
            actual_columns.append(column)
    print("... healthy samples removed")
    print(np.array([actual_columns]))
else:
    actual_columns = data.columns
data = data[actual_columns]
#data["accession"] = data_protein_names.iloc[:,0]
#data_clinical = data_clinical[data_clinical["ID_set"].isin(actual_columns)]
data_clinical = data_clinical.loc[data.columns, :]

## reverse log2
if log2 == True:
    data = 2. ** data

data["accession"] = ""
data["sort"] = np.nan
sample_columns = ~data.columns.isin(['sort', 'accession'])

## fill NA, prepare sorting, transform to zscore
for i, row in data.iterrows():
    mean = np.nanmean(data.iloc[i, sample_columns])
    SD = np.nanstd(data.iloc[i, sample_columns], ddof=0)
    if by_longlist == True:
        data.at[i, "sort"] = data_longlist.loc[data_longlist.iloc[:,0] == data_protein_names.iloc[i, 0]].iloc[0, 1]
        data.at[i, "accession"] = data_longlist.loc[data_longlist.iloc[:,0] == data_protein_names.iloc[i, 0]].iloc[0, 0]
    else:
        data.at[i, "sort"] = SD
        data.at[i, "accession"] = data_protein_names.iloc[i, 0]
    if zscore == True:
        data.iloc[i, sample_columns] = (data.iloc[i, sample_columns] - mean)/float(SD)
        data.iloc[i, sample_columns] = data.iloc[i, sample_columns].fillna(0)
    else:
        data.iloc[i, sample_columns] = data.iloc[i, sample_columns].fillna(mean)

### sort
data = data.sort_values(by="sort", ascending=False)
data = data.reset_index(drop=True)
data = data.iloc[int(num_proteins * cluster_quantile_upper):int(num_proteins * cluster_quantile_lower),:]
data.to_excel(data_folder + "/CLUSTER/_files/*input_data.xlsx", index=False)
print("... initial # of proteins:%s"%num_proteins)
print("... # proteins submitted to NMF (from %s to %s percentiles):"%(cluster_quantile_upper, cluster_quantile_lower),\
      len(data.iloc[:,0]),"(-%s %%)"%np.round((1-len(data.iloc[:,0])/float(num_proteins))*100,1))

### RUN NMF
np.save(data_folder + "/CLUSTER/_files/PP_patients",np.array(actual_columns))
data = np.array(data[actual_columns])
# level of sparsity
sparsity = (np.sum(data != 0)/data.size)
print("... sparsity:", (np.sum(data != 0)/data.size))
# negative elements
print("... negative elements", (np.sum(data < 0)), data.size)
# remove negative elements
data = data - np.min(data)
print("... negative elements after removal", (np.sum(data < 0)))
# normalize data
data = preprocessing.Normalizer().fit_transform(data)

def clean_axis(ax, y=False):
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().grid(False)
    if y == False: ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
            
fig = plt.figure(0, figsize=(13.9*.75*3*0.75, 10*.75*2.8*0.75))
outer = gridspec.GridSpec(3, 3, figure = fig, wspace=0.1, hspace=0.1, left=0.02)

cols = ["sparsity", "rank", "dispersion", "cophenetic"]
descriptives = pd.DataFrame(columns = cols)

print(">>> Starting NMF ...")
for i,rank in enumerate(range(2,11)):
    inner = gridspec.GridSpecFromSubplotSpec(1, 2,
                    subplot_spec=outer[i], wspace=0., hspace=0., width_ratios=[0.25,1])
    #######
    nmf = nimfa.Nmf(data, rank=rank,  max_iter=rank*100,
                    update='divergence', objective='div', n_run=100, seed='random_vcol', track_factor=True)
    nmf_fit = nmf()
    #######
    dispersion = nmf_fit.fit.dispersion()
    print("... dispersion k = ",rank,"\t", dispersion)
    C = 1 - nmf_fit.fit.consensus()
    #CC = scipy.spatial.distance.pdist(C)
    Y = sch.linkage(squareform(C), method='average')
    clusters = nmf_fit.fit.predict(prob=True) #what='features',
    clusters = np.array(clusters)

    # save data
    library = data_folder + "/CLUSTER/_files"

    np.save(library + "/PP_%s_k%s_C" % (mode, rank), C)

    np.save(library+"/PP_%s_k%s_clusters"%(mode,rank),clusters)
    clusters = pd.DataFrame(clusters)
    clusters.to_excel(library + "/PP_%s_k%s_clusters.xlsx" % (mode, rank), index=False)

    coefficients = np.array(nmf_fit.fit.coef())
    np.save(library+"/PP_%s_k%s_coefs"%(mode,rank),coefficients)

    coph = np.array(nmf_fit.fit.coph_cor())
    np.save(library+"/PP_%s_k%s_coph"%(mode,rank),coph)
    
    basis = np.array(nmf_fit.fit.basis())
    basis = pd.DataFrame(basis)
    basis.to_excel(library + "/PP_%s_k%s_basis.xlsx" % (mode, rank), index=False)

    scores = np.array(nmf_fit.fit.score_features())
    scores = pd.DataFrame(scores)
    scores.to_excel(library+"/PP_%s_k%s_scores.xlsx"%(mode,rank), index=False)

    selects = np.array(nmf_fit.fit.select_features())
    selects = pd.DataFrame(selects)
    selects.to_excel(library+"/PP_%s_k%s_selects.xlsx"%(mode,rank), index=False)

    descriptives = descriptives.append(pd.DataFrame([[sparsity, rank, dispersion, coph]], columns=cols))
    descriptives.to_excel(data_folder + "/CLUSTER/figures/descriptives.xlsx")

    # PLOT
    # combined figure
    denAX = fig.add_subplot(inner[0])
    denD = sch.dendrogram(Y, orientation='left', link_color_func=lambda k: 'black', no_labels=True)  # labels = labels
    clean_axis(denAX, y = True)
    heatmapAX = fig.add_subplot(inner[1])
    heatmapAX.set_title("k=%s"%rank)
    D = C[denD['leaves'], :][:, denD['leaves']]
    denD_temp = dict(denD)
    denD_temp = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in denD_temp.items()]))
    denD_temp.to_excel(library+"/PP_%s_k%s_denD.xlsx"%(mode,rank), index=False)
    axi = heatmapAX.imshow(D, interpolation='nearest', aspect='equal', origin='lower', cmap='RdBu') 
    clean_axis(heatmapAX)
    cb = fig.colorbar(axi, fraction=0.046, pad=0.04, aspect=10) 
    cb.set_label('Distance')#, fontsize=20)
    # rank-specific figure
    fig_single = plt.figure(rank, figsize=(13.9, 12)) #with labels 10 instead of 12.5
    heatmapGS = gridspec.GridSpec(len(cohort_features)+1, 2, figure = fig_single, wspace=.0, hspace=0.018, width_ratios=[0.25, 1], height_ratios=height_ratios)
    denAX = fig_single.add_subplot(heatmapGS[len(cohort_features), 0])
    denD = sch.dendrogram(Y, orientation='left', link_color_func=lambda k: 'black', no_labels=True)
    clean_axis(denAX, y=True)
    heatmapAX = fig_single.add_subplot(heatmapGS[len(cohort_features), 1])
    axi = heatmapAX.imshow(D, interpolation='nearest', aspect='auto', origin='lower', cmap='RdBu')
    clean_axis(heatmapAX)
    # features
    for j,l in enumerate(cohort_features):
        featureAX = fig_single.add_subplot(heatmapGS[j, 1])
        bar = featureAX.imshow([data_clinical.loc[actual_columns, l].values[denD['leaves']]], cmap=cmaps[j], \
                               vmin = cmap_low[j], vmax= cmap_up[j], aspect='auto')
        featureAX.text(-0.01, .5, l, transform=featureAX.transAxes,
                size=15, ha='right', va='center') #weight='bold')
        clean_axis(featureAX)
    fig_single.savefig(data_folder + "/CLUSTER/figures/PP_%s_clusters_k%s.png" % (mode,rank), dpi=300)
fig.savefig(data_folder + "/CLUSTER/figures/PP_%s_clusters_all.png"%(mode),dpi=300)