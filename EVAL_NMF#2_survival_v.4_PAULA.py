import numpy as np
import pylab as plt
import sys, os, inspect
directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append('/Users/dressler/Documents/Core/_PythonModules')
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
#print data_survival

###
cls, cls_prob = np.load(data_folder + "/CLUSTER/_files/PP_%s_k%s_clusters.npy"%(mode,k)) #needs pip install numpy==1.16.1 for this to work
cls = np.array(cls)[0]
print len(cls)
pts = np.load(data_folder + "/CLUSTER/_files/PP_patients.npy")
pts = list(pts)
if len(cls) != len(pts):
    pts = pts[1:]
if len(cls) != len(pts):
    print "ERROR: Check patients and cluster lists - incompatible lengths"
    sys.exit()

###
if by_feature != None:
    cls = []
    for pt in pts:
        cls.append(int(data_clinical.where(data_clinical["ID_set"] == pt).loc[:,by_feature].dropna().iloc[0]))
    cls = np.array(cls)

###


data_cls = []
# [time observed, event (1,0)
for l in range(len(np.unique(cls))): data_cls.append([[],[]])
fail = 0
for i,c in enumerate(cls):
    if c == combine[0]:
        c = combine[1]
    pt = pts[i]
    found = False
    for j,p in enumerate(data_survival.loc[:,id_column]):
        if "U" in p and p in pt and exclude not in pt:
            print "found",pt
            data_cls[c][0].append(data_survival.loc[j, duration_column])
            data_cls[c][1].append(data_survival.loc[j, event_column])
            found = True
    if found == False:
        print("... %s not found"%pt)
        fail+=1
print "not found", fail
    
###
fig = plt.figure(1,[8,7])
ax = fig.add_subplot(111)
Sum = 0
os = []
for i,d_cls in enumerate(data_cls):
    os.append([])
    number = len(d_cls[0])
    Sum += number
    print(number)
    durations, event_observed = d_cls
    for d,e in zip(durations, event_observed):
        if e == 1:
            os[i].append(d)
    km = KaplanMeierFitter()
    try:
        km.fit(durations, event_observed, label='Cluster #%s'%(i))
        km.plot(ax = ax, ci_show = False, show_censors = True) #, color = colors[i])
    except:
        print durations, event_observed
print("N=%s"%Sum)


import statistics as ss
for i,s in enumerate(os):
    print ""
    print ""
    try:
        ss.sumup(s, name="cluster #%s" % (i))
    except:
        pass
    for j in range(i+1,len(os)):
        print ""
        print "... cluster #%s vs. cluster #%s" %(i,j)
        try:
            #ss.deltaquant(os[i],os[j],names = ["cluster #%s" %(i),"cluster #%s" %(j)])
            lstat.logrank_test(data_cls[i][0], data_cls[j][0], data_cls[i][1], data_cls[j][1]).print_summary()
        except: pass
#ax.set_xlim(0,100)
ax.set_ylim(0,1)
ax.set_xlabel("survival [months]")
ax.set_ylabel("share of patients at risk [1]")
ax.grid(False)

if by_feature != None:
    plt.savefig(data_folder + "/CLUSTER/figures/survival_by-%s_k%s.png" %(by_feature,k), dpi=300)
else:
    plt.savefig(data_folder + "/CLUSTER/figures/survival_k%s.png"%k,dpi=300)
plt.show()