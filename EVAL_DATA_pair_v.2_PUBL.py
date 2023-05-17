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
vs.style('science')
colors = vs.colors_indexed

plot = bool(raw_input("... plot individual distributions? [True/False] "))

def p_value(fold_sigma):
    """one-sided; modified one-sample ttest"""
    return 1. - stats.norm.cdf(fold_sigma)

def create_null(primary,all):
    """name of primary healthy sample, list of all healthy samples"""
    global data_temp
    temp = data_temp[all]
    temp = temp.where(temp>-10,other=np.nan)
    for secondary in all:
        if secondary != primary:
            temp[secondary] -= temp[primary]
    temp.drop(primary,'columns',inplace=True)
    return temp


### load normalized data
data_folder = "*PAULA"
data = pd.read_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICSmed-log2_COMPLETE-DATA.xlsx")
data["Accession"] = pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins_COMPLETE-DATA.xlsx")["Accession"]
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical.xlsx")
data_functional = pd.read_excel(data_folder + "/*RES_data_functional.xlsx")

### set samples
samples = []
for sample in data.columns:
    if "U-" in sample:
        samples.append(sample)

### iterate over TMT sets
columns_functional = ["# unique PSM"]
data_paired = pd.DataFrame()
data_paired["Accession"] = data["Accession"]
data_significance = pd.DataFrame()
data_significance["Accession"] = data["Accession"]
try:
    os.mkdir(data_folder + "/PAIR")
except:
    pass
if plot == True:
    try:
        os.mkdir(data_folder + "/PAIR/figures_individual")
    except:
        pass

for TMT_set in data_clinical["TMT_set"].unique(): #["H01"]
    print("... pairing samples from set #%s"%TMT_set)
    ## define TMT set subset
    samples_set = []
    samples_set_healthy = []
    samples_set_tumor = []
    for sample in samples:
        if list(data_clinical.where(data_clinical["ID_set"] == sample).loc[:,"TMT_set"].dropna())[0] == TMT_set:
            samples_set.append(sample)
            if "_H" in sample:
                samples_set_healthy.append(sample)
            else:
                samples_set_tumor.append(sample)
    ## collect peptide data
    columns_functional_set_specific = []
    for column in columns_functional:
        columns_functional_set_specific.append(TMT_set + "_" + column)
    data_temp = pd.merge(data[["Accession"] + samples_set], data_functional[["Accession"] + columns_functional_set_specific], how='left', on=["Accession"])

    ## iterate over tumor samples
    for sample_tumor in samples_set_tumor:
        # ensure pairs
        if sample_tumor[:-1]+"H" in samples_set_healthy:
            # load samples for null distribution
            data_null = create_null(sample_tumor[:-2]+"_H",samples_set_healthy)
            # create peptide score based on #_unique_peptides ** 2 + #_of_other_peptides
            data_null["Peptide confidence score"] = (data_temp[TMT_set + "_" + "# unique PSM"])
            # normalize healthy to tumor sample
            data_null[sample_tumor + "_paired"] = data_temp[sample_tumor] - data_temp[sample_tumor[:-2]+"_H"]
            # deal with pairs where either value was zero (were replaced with 10**-6 during normalization because of log transformation, create log2 ratios above or below 10)
            min = np.nanmin(data_null[sample_tumor + "_paired"].where(data_null[sample_tumor + "_paired"] > -10))
            max = np.nanmax(data_null[sample_tumor + "_paired"].where(data_null[sample_tumor + "_paired"] < 10))
            paired_values = []
            for T,H in zip(data_temp[sample_tumor], data_temp[sample_tumor[:-2]+"_H"]):
                if math.isnan(T) == False and math.isnan(H) == False:
                    if T-H < -10.:
                        paired_values.append(min)
                    elif T-H > 10.:
                        paired_values.append(max)
                    else:
                        paired_values.append(T-H)
                else:
                    paired_values.append(np.nan)
            data_null[sample_tumor + "_paired"] = paired_values
            # PLOT tumor versus null distribution across peptide score range
            if plot == True:
                plt.close('all')
                fig = plt.figure("Check %s"%sample_tumor, [12,12])
                fig.subplots_adjust(hspace=.4, wspace=.3)
                ax = fig.add_subplot(4,4,4)
                ax.grid(False)
                ax.set_xlabel("Unique PSMs [log$_2$]")
                ax.set_ylabel("Fold change [log$_2$]")
                ax.plot(np.log2(data_null["Peptide confidence score"]), data_null[sample_tumor + "_paired"], ".",
                        alpha=.1, color = colors[1],markeredgewidth=.0, label="_nolegend_")
                for sample in samples_set_healthy:
                    if sample_tumor[:-1] not in sample:
                        #Z, xedges, yedges = np.histogram2d(np.log(data_null["Peptide confidence score"]), data_null[sample])
                        #ax.pcolormesh(xedges, yedges, Z.T)
                        ax.plot(np.log2(data_null["Peptide confidence score"]), data_null[sample], ".",
                                label="_nolegend_", alpha=.1/6., color=colors[0], markeredgewidth=.0)
            # calculate x (= peptide score) specific variances to define local null distribution
            xs, ys, bl, bh, p_low, p_up = [],[],[],[], [], []
            for x in data_null["Peptide confidence score"].sort_values(ascending=True).unique():
                # bins are of cumulative size to compensate for decreasing density of values
                bin_low = .5 * x
                bin_high = 2. *  x
                bl.append(bin_low)
                bh.append(bin_high)
                bin = data_null.where((data_null["Peptide confidence score"] >= bin_low) & (data_null["Peptide confidence score"] <= bin_high)).loc[:, data_null.columns != "Peptide confidence score"]
                bin = bin.to_numpy().flatten()
                xs.append(x)
                ys.append(np.nanvar(bin)**.5)
                p_low.append(np.nanpercentile(bin, 2.5))
                p_up.append(np.nanpercentile(bin, 97.5))
            # PLOT upper and lower significance thresholds
            if plot == True:
                #ax.plot(np.log2(xs), ys, "--", color="grey")
                ax.plot(np.log2(xs), np.array(ys) * 1.65, "--", color="black")
                ax.plot(np.log2(xs), np.array(ys) * -1.65, "--", color="black")
                #ax.plot(np.log2(xs), np.array(ys) * -1., "--", color="grey")
                #ax.plot(np.log2(xs), np.array(p_low), "--", color=colors[1])
                #ax.plot(np.log2(xs), np.array(p_up), "--", color=colors[1])
                ax.plot([None],[1.],".", marker="s", color=colors[0], label="null distr.", alpha=.75)
                ax.plot([None],[1.],".", marker="s", color=colors[1], label="sample", alpha=.75)
                ax.set_title(sample_tumor+"\n")
                ax.legend(fancybox=False)
                fig.savefig(data_folder + "/PAIR/figures_individual/%s.png" % sample_tumor, dpi=300, transparent=True)
            library = {}
            for x,y in zip(xs,ys):
                library[x] = y
            # calculate p- and q-values
            data_p_values = []
            for i,row in data_null.iterrows():
                x = row["Peptide confidence score"]
                y = row[sample_tumor + "_paired"]
                if math.isnan(x) == False and math.isnan(y) == False:
                    pval = p_value(np.abs(y)/float(library[x]))
                    data_p_values.append(pval)
                else:
                    data_p_values.append(np.nan)
            data_p_values = np.array(data_p_values)
            q_values = multitest.fdrcorrection(data_p_values[~np.isnan(data_p_values)], alpha=0.05)[1]
            # refill q_values with nan values to match other columns
            data_q_values = []
            q_counter = 0
            for p in data_p_values:
                if math.isnan(p) == False:
                    data_q_values.append(q_values[q_counter])
                    q_counter += 1
                else:
                    data_q_values.append(np.nan)
            print len(q_values), q_counter
            # append to data
            data_paired[sample_tumor + "_paired"] = data_null[sample_tumor + "_paired"]
            data_significance[sample_tumor + "_p-value"] = data_p_values
            data_significance[sample_tumor + "_q-value"] = data_q_values


    #fig = plt.figure("Binning function")
    #ax = fig.add_subplot(111)
    #ax.grid(False)
    #ax.plot(np.log(xs),np.log(bl),"--",color=colors[0])
    #ax.plot(np.log(xs),np.log(xs),color="black")
    #ax.plot(np.log(xs),np.log(bh),"--",color=colors[0])
    #plt.show()

data_paired.to_excel(data_folder + "/PAIR/*RES_data_paired.xlsx",index=False)
data_significance.to_excel(data_folder + "/PAIR/*RES_data_significance.xlsx",index=False)