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
import statistics as st
import string
import random
from sklearn.decomposition import PCA
vs.style('science')
colors = vs.colors_indexed

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)


### settings
data_folder = "*PAULA"
select = 5

data = pd.read_excel(data_folder + "/CLUSTER/figures/descriptives.xlsx")

fig = plt.figure("Rank determination", [10, 10])
fig.subplots_adjust(hspace=.4, wspace=.4)

ax = fig.add_subplot(332)
#ax.set_title("Dispersion")
ax.set_ylabel("Dispersion [1]")
ax.set_xlabel("Number of clusters (k) [1]")
ys = data.loc[:,"dispersion"]
ax.plot([select, select],[min(ys),max(ys)],"--", color="black")
ax.plot(data.loc[:,"rank"],ys, color=colors[0])
ax.grid(False)

ax = fig.add_subplot(331)
#ax.set_title("Cophenetic correlation")
ax.set_ylabel("Cophenetic correlation [1]")
ax.set_xlabel("Number of clusters (k) [1]")
ys = data.loc[:,"cophenetic"]
ax.plot([select, select],[min(ys),max(ys)],"--", color="black")
ax.plot(data.loc[:,"rank"],ys, color=colors[0])
ax.grid(False)

plt.show()
fig.savefig(data_folder + "/CLUSTER/figures/PLOT_rank-determination.png", dpi=300)
