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
#min_up = 10
max_up = 150
min_up = 10
max_dn = 150
min_dn = 10
sig = .05
sig_mode = "p" # either "p" (FDR off) or "q" (FDR on)
min_correlation = .5
sig_correlation = .05
min_FC_up = .25
min_FC_dn = -.25
#########################

try: os.mkdir(data_folder + "/CMAP")
except: pass
try: os.mkdir(data_folder + "/CMAP/sets")
except: pass

### retrieve current l1000 gene list
genes_available = list(pd.read_excel(data_folder + "/CMAP/data_l1000-genes_BING.xlsx", dtype='string')["gene_id"])

### retrieve available mRNA-protein correlations
correlations_available = pd.read_excel(data_folder + "/CMAP/data_correlation.xlsx").loc[:,["Gene name","Spearman correlation","p-value"]]

### load pairwise normalized data
data = pd.merge(\
    pd.read_excel(data_folder + "/PAIR/*RES_data_paired.xlsx"), \
    pd.read_excel(data_folder + "/*RES_data_functional.xlsx"), \
    how='left', on=["Accession"])
data_significance = pd.read_excel(data_folder + "/PAIR/*RES_data_significance.xlsx")

samples = []
for sample in data.columns:
    if "U-" in sample:
        samples.append(sample)

### loop through samples
res = pd.DataFrame()
for i in range(len(samples)):  #range(2): 
    sample = samples[i]

    ## sort data
    temp = data.loc[:,["Entrez Gene ID","Gene Symbol","Accession",sample]]
    temp = temp.dropna()
    column_sig = "%s_T_%s-value"%(sample.split("_")[0],sig_mode)
    temp = pd.merge(temp, data_significance[["Accession",column_sig]], how='left', on=["Accession"])
    temp = temp.sort_values(by=column_sig, ascending=True)

    ## find up-regulated set
    up = []
    up.append(sample) # name
    up.append("") # optional description of the sample (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
    for u,row in temp.iterrows():
        #if len(up)-2 >= min_up \
         #       and len(up)-2 < max_up \
          #      and row[sample] > min_FC_up \
           #     and row["Entrez Gene ID"] in genes_available \
            #    and row[column_sig] < sig:
            #try:
            #    correlation = correlations_available.where(correlations_available["Gene name"] == row["Gene Symbol"]).dropna().iloc[0,:]
            #    r = correlation["Spearman correlation"]
            #    p = correlation["p-value"]
            #    if r > min_correlation and p <= sig_correlation: # and len(are_genes_in_api(client,[row["Gene Symbol"]])) == 1:
            #        up.append(row["Entrez Gene ID"])
            #except:
            #    pass
        if len(up)-2 < max_up \
                and row[sample] > min_FC_up \
                and row["Entrez Gene ID"] in genes_available\
                and row[column_sig] < sig:
            try:
                correlation = correlations_available.where(correlations_available["Gene name"] == row["Gene Symbol"]).dropna().iloc[0,:]
                r = correlation["Spearman correlation"]
                p = correlation["p-value"]
                if r > min_correlation and p <= sig_correlation: # and len(are_genes_in_api(client,[row["Gene Symbol"]])) == 1:
                    up.append(row["Entrez Gene ID"])
            except:
                pass
    print("... sample %s: UP with %s length"%(sample, len(up)))

    ## find down-regulated set
    #temp = temp.sort_values(by=sample, ascending=True)
    dn = []
    dn.append(sample) # name
    dn.append("") # optional description https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
    for u,row in temp.iterrows():
        if len(dn)-2 < max_dn \
                and row[sample] < min_FC_dn\
                and row["Entrez Gene ID"] in genes_available\
                and row[column_sig] < sig: # and row["Entrez Gene ID"] in genes_available:
            try:
                correlation = correlations_available.where(correlations_available["Gene name"] == row["Gene Symbol"]).dropna().iloc[0,:]
                r = correlation["Spearman correlation"]
                p = correlation["p-value"]
                if r > min_correlation and p <= sig_correlation: # and len(are_genes_in_api(client,[row["Gene Symbol"]])) == 1:
                    dn.append(row["Entrez Gene ID"])
            except:
                pass
    print("... sample %s: DN with %s length"%(sample, len(dn)))

    k = i/25
    if len(up)-2 >= min_up and len(dn)-2 >= min_dn:
        incl = 1
        with open(data_folder + "/CMAP/sets/batch_%s_UP.gmt"%k,"a") as f:
            f.write("\t".join(up))
            if np.mod(i+1,25) != 0: f.write("\n")
        with open(data_folder + "/CMAP/sets/batch_%s_DN.gmt"%k,"a") as f:
            f.write("\t".join(dn))
            if np.mod(i + 1, 25) != 0: f.write("\n")
    else:
        incl = 0
        print("... skipped sample %s"%sample)
    res = res.append(pd.DataFrame([[sample, len(up), len(dn), incl, k]], columns=["sample","UP", "DN", "Included", "batch"]))
res.to_excel(data_folder + "/CMAP/INFO_batch-summary.xlsx", index=False)