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
cluster = int(raw_input("... cluster k= "))
cluster_order = str(raw_input("... order [int separated by ,] = "))
exec("cluster_order = [%s]"%cluster_order)
print("... selected order was: %s"%cluster_order)
mode_index = 2
extract_image_order = True
#########################
library = data_folder + "/CLUSTER/_files"
modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedNORMICSmed-log2_PM",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
mode = modes[mode_index]
clusters = np.array(np.load(library+"/PP_%s_k%s_clusters.npy"%(mode,cluster))[0])[0]
samples = np.load(library+"/PP_patients.npy")
print clusters
print samples

data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical.xlsx", index_col="ID_set")

if extract_image_order == True:
    denD = pd.read_excel(library+"/PP_%s_k%s_denD.xlsx"%(mode,cluster))
    denD = denD.sort_values(by="leaves", ascending=True)
    image_order = list(denD.index)
    for o, s in zip(image_order, samples):
        print("... extracting consensus plot order %s for sample %s"%(o,s))
        data_clinical.at[s, "image_order k=%s" % cluster] = o

for c,s in zip(clusters,samples):
    print("... extracting cluster %s for sample %s"%(c,s))
    data_clinical.at[s,"cluster k=%s"%cluster] = c
data_clinical["cluster k=%s"%cluster] = data_clinical["cluster k=%s"%cluster].replace(\
    cluster_order, np.sort(data_clinical["cluster k=%s"%cluster].dropna().unique()))
data_clinical.to_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")
