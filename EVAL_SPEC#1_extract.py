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
vs.style('science')
colors = vs.colors_indexed
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

#######################
data_folder = "*PAULA"
sig = 0.05
sig_mode = "_p-value"
#######################

data = pd.read_excel(data_folder + "/PAIR/*RES_data_paired.xlsx").iloc[:,:]
data = pd.merge(data, pd.read_excel(data_folder + "/*RES_data_functional.xlsx"), how='left', on='Accession')
data_sig = pd.read_excel(data_folder + "/PAIR/*RES_data_significance.xlsx").iloc[:,:]
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")

res = pd.DataFrame()
res["Accession"] = data["Accession"]
res["Gene Symbol"] = data["Gene Symbol"]
res[["Description","Biological Process","Entrez Gene ID","Ensembl Gene ID"]] = data[["Description","Biological Process","Entrez Gene ID","Ensembl Gene ID"]]
groups_stage = ["ALL","pTa low", "pTa high", "pTis", "pT1a/b/G1/G2", "pT1c/G3", "pT2", ">pT2","pN1"]
groups_cluster = ["PAULA I", "PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III"]
conditions = []
for group, condition in zip(groups_stage, [">0"] + ["==%s"%x for x in range(1,9)]):
    conditions.append((group, "Stage",condition))
for group, condition in zip(groups_cluster, ["==%s"%x for x in range(0,5)]):
    conditions.append((group,"cluster k=5",condition))

samples_ori = [x for x in data.columns if x[:2] == "U-"]

print(">>> Start specificity extraction ...")
for Set in conditions:
    group, parameter, condition = Set
    print("... calculating for %s"%group)

    res["%s_direction" % (group)] = np.nan
    for direction in ["UP", "DN"]:
        res["%s_%s_n"%(group,direction)] = np.nan
        res["%s_%s_share"%(group,direction)] = np.nan
        res["%s_%s_share-available"%(group,direction)] = np.nan
        res["%s_%s_median"%(group,direction)] = np.nan
        res["%s_%s_mean"%(group,direction)] = np.nan
        res["%s_%s_min"%(group,direction)] = np.nan
        res["%s_%s_max"%(group,direction)] = np.nan
        res["%s_%s_p25"%(group,direction)] = np.nan
        res["%s_%s_p75"%(group,direction)] = np.nan

    exec('samples_possible_ori = data_clinical["ID_set"].where(data_clinical["%s"]%s).dropna().values'%(parameter,condition))

    samples_possible = [x for x in samples_ori if x[:-7] in samples_possible_ori]

    for i,D in data.iterrows():
        samples = np.array(samples_possible)[~D[samples_possible].isna()]
        samples_p = [x[:-7] + sig_mode for x in samples]
        for direction in ["UP","DN"]:
            if direction == "UP":
                samples_sig = [x <= sig and y > 0 for x,y in zip(data_sig.loc[i, samples_p].values,
                                                       D[samples].values)]
            else:
                samples_sig = [x <= sig and y < 0 for x, y in zip(data_sig.loc[i, samples_p].values,
                                                                  D[samples].values)]
            samples_count = len([x for x in samples_sig if x != False])
            samples_values = D[samples].values[samples_sig]
            if samples_count > 0:
                res.at[i,"%s_%s_n" % (group, direction)] = samples_count
                res.at[i,"%s_%s_share" % (group, direction)] = samples_count/float(len(samples_possible))
                res.at[i,"%s_%s_share-available" % (group, direction)] = samples_count/float(len(samples))
                res.at[i,"%s_%s_median" % (group, direction)] = np.median(samples_values)
                res.at[i,"%s_%s_mean" % (group, direction)] = np.mean(samples_values)
                res.at[i,"%s_%s_min" % (group, direction)] = np.min(samples_values)
                res.at[i,"%s_%s_max" % (group, direction)] = np.max(samples_values)
                res.at[i,"%s_%s_p25" % (group, direction)] = np.percentile(samples_values, 25.)
                res.at[i,"%s_%s_p75" % (group, direction)] = np.percentile(samples_values, 75.)
            else:
                res.at[i, "%s_%s_n" % (group, direction)] = 0
                res.at[i, "%s_%s_share" % (group, direction)] = 0.
                res.at[i, "%s_%s_share-available" % (group, direction)] = 0.
        res.at[i,"%s_direction" % (group)] = np.log2((res.at[i,"%s_UP_n" % (group)]+1)
                                                     / (float(res.at[i,"%s_DN_n" % (group)])+1.))

try:
    os.mkdir(data_folder + "/SPECIFICITY")
except:
    pass
res.to_excel(data_folder + "/SPECIFICITY/*RES_specificity.xlsx", index=False)