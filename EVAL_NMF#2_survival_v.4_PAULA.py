import numpy as np
import pylab as plt
import sys, os, inspect
directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
sys.path.append(file_directory + '/_modules')
import visualize as vs
import xml.etree.ElementTree as elt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines import statistics as lstat

vs.style('science')
colors = vs.colors_indexed

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)



### set basic parameters ###
data_folder = "*PAULA"
k = raw_input("k = ")
mode = 2
id_column = "ID"
duration_column = "OS"
event_column = "Event"
exclude = "_H" #"?" #"_H"
combine = [None, None]#None,None]
by_feature = None #None
############################



modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedNORMICSmed-log2_PM",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
mode = modes[mode]
#data_z = pd.read_excel(data_folder + "/CLUSTER/RES_data_%s_zscore-nolog.xlsx"%mode[:-5])
#data = pd.read_excel(data_folder + "/CLUSTER/RES_data_%s_nolog.xlsx"%mode[:-5])
if "_PM" in mode:
    data_protein_names = pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins_PM.xlsx")
else:
    data_protein_names = pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins.xlsx")
num_proteins = len(data_protein_names.iloc[:, 0])
data_survival = pd.read_excel(data_folder + "/*PAULA_survival.xlsx")
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical.xlsx")
