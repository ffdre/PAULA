### import modules
import urllib
import numpy as np
import pylab as plt
import pandas as pd
import multiprocessing as mp
import subprocess
import sys, os, inspect
import scipy
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
sys.path.append(file_directory + '/_modules')
import visualize as vs
import statistics as st
import random
import string
from operator import itemgetter

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

################################
normalize_batch = True
fuse_first = True
################################

if fuse_first == True:
    data = pd.read_excel("*RES_data_raw_normalize_batch=%s.xlsx"%normalize_batch)
else:
    data = pd.read_excel("*RES_data_raw_normalize_batch=%s_fused.xlsx" % normalize_batch)
instructions = pd.read_excel("*columns.xlsx")
columns_to_fuse = list(instructions["columns_to_fuse_in_clean-up"].dropna())
sample_columns = []
for column in data.columns:
    if column[:2] == "U-":
        sample_columns.append(column)

### fuse columns
if fuse_first == True:
    for column_to_fuse in columns_to_fuse:
        print("... fusing %s"%column_to_fuse)
        set_to_fuse = []
        for column in data.columns:
            if column_to_fuse == column[4:]:
                set_to_fuse.append(column)
                print(column)
        new_column = []
        for i,row in data[set_to_fuse].iterrows():
            row_fused = list(row.dropna())
            if len(row_fused) > 0:
                new_column.append(row_fused[0])
            else:
                new_column.append(None)
        # check column lenght
        if len(new_column) != len(data.iloc[:,0]):
            raise ValueError
        data.drop(set_to_fuse, "columns", inplace=True)
        data[column_to_fuse] = new_column
    data.to_excel("*RES_data_raw_normalize_batch=%s_fused.xlsx"%normalize_batch, index=False)

### remove proteins with no intensities
print("... %s proteins originally"%len(data.iloc[:,0]))
data.dropna(subset=sample_columns, inplace=True, thresh=1)
data.reset_index(drop=True, inplace=True)
print("... %s proteins after removal of proteins without any intensities"%len(data.iloc[:,0]))

### remove isoforms
dict_of_excluded = {}
dict_of_excluded_accession = {}
excluded_all = []
checked = []

def check_for_isoform(i):
    global data, sample_columns, dict_of_excluded, dict_of_excluded_accession, initial_length, excluded_all, checked
    data_i = data.loc[i,sample_columns]
    multiples, multiples_names, multiples_ids = [], [], []
    for j in range(i,initial_length):
        if data.loc[j,sample_columns].equals(data_i) == True:
            if j == i+1:
                print("... detected duplicate in %s" % data.loc[i, "Description"])
            if j > i:
                print("... ... duplicate         %s" % data.loc[j, "Description"])
                checked.append(j)
            multiples.append(j)
            multiples_names.append(data.loc[j, "Description"])
            multiples_ids.append(data.loc[j, "Accession"])
        else:
            break
    if len(multiples) > 1:
        #print multiples_names
        #keep = multiples[int(raw_input("... choose index of isoform to keep: "))]
        excluded, excluded_names, excluded_accessions = [], [], []
        keep = 0
        for k,m in enumerate(multiples):
            if multiples_ids[k][0] == "P":
                keep = k
            elif multiples_ids[k][0] != "A":
                keep = k
        for k, m in enumerate(multiples):
            if k != keep:
                excluded_all.append(m)
                excluded.append(m)
                excluded_names.append(multiples_names[k])
                excluded_accessions.append(multiples_ids[k])

        dict_of_excluded[multiples_names[keep]] = pd.Series(excluded_names)
        dict_of_excluded_accession[multiples_ids[keep]] = pd.Series(excluded_accessions)
        pass
    else:
        pass

data.sort_values(by=sample_columns[0], inplace=True)
initial_length = len(data.iloc[:,0])
for i in range(initial_length):
    if i not in checked:
        check_for_isoform(i)
print("... dropping duplicates")
data.drop(list(set(excluded_all)), "index", inplace=True)
dict_of_excluded = pd.DataFrame(dict_of_excluded)
dict_of_excluded.to_excel("*INFO_data_raw_normalize_batch=%s_clean-excluded.xlsx"%normalize_batch, index=False)
dict_of_excluded_accession = pd.DataFrame(dict_of_excluded_accession)
dict_of_excluded_accession.to_excel("*INFO_data_raw_normalize_batch=%s_clean-excluded-accessions.xlsx"%normalize_batch, index=False)
print("... saving cleaned-up data")
#data.to_csv("*RES_data_raw_normalize_batch=%s_clean.csv"%normalize_batch, index=False)
data.to_excel("*RES_data_raw_normalize_batch=%s_clean.xlsx"%normalize_batch, index=False)

print("... all done and saved. <<<")

