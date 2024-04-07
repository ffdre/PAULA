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
from adjustText import adjust_text
import statsmodels.stats.api as sms
vs.style('science')
colors = vs.colors_indexed

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)


### settings
data_folder = "*PAULA"
modes = ["non-normalized-log2",
         "normalizedNORMICS-glog",
         "normalizedNORMICSmed-log2",
         "normalizedVSN-glog",
         "normalizedMEDIAN-log2",
         "normalizedQUANTILE-log2",
         "normalizedLOESS-log2"]
titles = ["Non-normalized",
          "Normics",
          "Normics$_{median}$",
          "VSN",
          "Median",
          "Quantile",
          "Cyclic LOESS",
          "Significant DE proteins",
          "Fold-changes (DE)"]
mode_index = 2
plot_ci = True #log2_for_test must be True then!
FDR = True
MWU = True
log2_for_test = True # not relevant for MWU only for ttest
code = None #"ECOLI"
limit_x = 0.5 #None#np.log2(1.25)#.25#None #.5#.75#None #np.log2(1.5)
limit_y = -np.log10(0.05)
minimum_values = 0.1
group_column = "cluster k=5"
ID_column = "ID_set"
#Group1 = 1
Group2 = -1
group_colors = [colors[0], colors[1]] # different colors for group1, group2 significance
subplots_x = 2
subplots_y = 2
annotations = False
fontsize_annotations = 7
separator = " "
max_names = 3
y_min = None
y_max = None
subplot_letter = ""
markersize_sign = 5 #12.  #standard 7., 10. for DE only
markersize = 3 #10.


if FDR == True: titles.append("q-Values (DE)")
else: titles.append("p-Values (DE)")
if group_colors != None:
    colors = ["grey", group_colors[0], group_colors[1]]
else:
    colors = [colors[0], colors[1], colors[1]]

# allocate patients to groups
protein_names = pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins.xlsx").iloc[:, 0] #_annotated
data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")
data_clinical["cluster k=5"] = data_clinical["cluster k=5"].replace(np.nan,-1)
data_file = data_folder + "/NORMICS_results/RES_data_%s.xlsx"%modes[0]
protein_names = pd.merge(pd.read_excel(data_folder + "/NORMICS_results/RES_data_proteins.xlsx"),
                         pd.read_excel(data_folder + "/*RES_data_functional.xlsx")[["Accession", "Gene Symbol"]],
                         how='left', on='Accession')
protein_names.to_excel(data_folder + "/NORMICS_results/CTRL_data_proteins.xlsx", index=False)
protein_names = protein_names["Gene Symbol"].values
print protein_names
data = pd.read_excel(data_file)

patients = data.columns[:]
groups = {} #dict containing patients per group
group_list = []
for group in np.sort(data_clinical.loc[:, group_column].unique()):
    groups[group] = []
    group_list.append(group)
print group_list

for patient in patients:
    for ID, group in zip(data_clinical.loc[:, ID_column], data_clinical.loc[:, group_column]):
        if ID == patient:
            groups[group].append(patient)
        elif data_folder == "TCGA" and ID in patient:
            groups[group].append(patient)


### plot
eval_num = []
eval_change = []
eval_qvals = []
eval_var = []
subplot = 0
mode = modes[mode_index]
for Group1 in range(0,5):
    plt.close()
    # fig = plt.figure("Volcano plot(s) %s" % data_folder, [8, 8]) #10, 10
    fig = plt.figure("PLOT", [12, 12])
    fig.subplots_adjust(hspace=.4, wspace=.3)
    labels_xs, labels_ys, labels = [], [], []
    eval_change.append([])
    eval_qvals.append([])
    eval_var.append([])
    ax = fig.add_subplot(4,4,subplot+1)
    ax.text(-0.1, 1.05, subplot_letter, transform=ax.transAxes,
            size=20, weight='bold')
    #ax.set_title(titles[mode_index])

    data_file = data_folder + "/NORMICS_results/RES_data_%s.xlsx"%mode
    data = pd.read_excel(data_file)
    num_proteins = len(data.iloc[:, 0])

    # calculate DE
    group_data = {} #dict containing data per group
    for group in groups.keys():
        group_data[group] = []
    group_data["protein"] = []

    for i,p in enumerate(protein_names):

        values = list(data.iloc[i,:])
        for group in groups.keys():
            group_data[group].append([])
            for j,patient in enumerate(patients):
                if patient in groups[group]:
                    group_data[group][-1].append(values[j])

    for group in groups.keys():
        print "... loaded group '%s' with %s values"%(group,len(group_data[group][0]))

    for j, group1 in enumerate([Group1]):
        for k, group2 in enumerate([Group2]):
            k += (j+1)
            # calculate goupwise log-fold change and p-value


            res = pd.DataFrame()
            x_col = "fold_change_plot-cutoff=%s"%limit_x
            y_col = "p-value_plot-cutoff=%s"%np.round(10.**(-limit_y),3)
            cols = ["protein", x_col, y_col, "N1", "N2"]



            above_sign = 0
            left, right = 0, 0
            norm = 0
            pvals, X, N1, N2, CI = [], [], [], [], []
            for i,p in enumerate(protein_names):
                A = np.array(group_data[group1][i])
                if log2_for_test == False:
                    A = 2. ** A
                A = [a for a in A if (math.isnan(a) == False)]
                B = np.array(group_data[group2][i])
                if log2_for_test == False:
                    B = 2. ** B
                B = [b for b in B if (math.isnan(b) == False)]
                N1.append(len(A))
                N2.append(len(B))
                if len(A) > (minimum_values * len(group_data[group1][i])) \
                        and len(B) > (minimum_values * len(group_data[group2][i])):
                    if MWU == False:
                        pval = stats.ttest_ind(A, B)[1]
                        if math.isnan(pval) == False:
                            pvals.append(pval)
                    else:
                        pval = stats.mannwhitneyu(A, B)[1]
                        if math.isnan(pval) == False:
                            pvals.append(pval)
                    if math.isnan(pval) == False:
                        A = np.array(group_data[group1][i])
                        A = [a for a in A if (math.isnan(a) == False)]
                        B = np.array(group_data[group2][i])
                        B = [b for b in B if (math.isnan(b) == False)]
                        X.append(np.mean(A) - np.mean(B))  # np.log2(np.mean(A) / float(np.mean(B)))
                        if plot_ci == True:
                            cm = sms.CompareMeans(sms.DescrStatsW(A), sms.DescrStatsW(B))
                            ci = cm.tconfint_diff(usevar='unequal')
                            CI.append(ci)
                    else:
                        X.append('NA')
                        if plot_ci == True:
                            CI.append('NA')
                else:
                    X.append('NA')
                    if plot_ci == True:
                        CI.append('NA')
            #print pvals
            if FDR == True:
                qvals = multitest.fdrcorrection(np.array(pvals), alpha=0.05)[1]
            else:
                qvals = pvals
            i_qvals = 0
            for i, p in enumerate(protein_names):
                x = X[i]
                ci = CI[i]
                n1, n2 = N1[i], N2[i]
                if x != 'NA':
                    y = -np.log10(qvals[i_qvals])
                    i_qvals += 1
                    if code == None:
                        #if "norm" in p:
                        #    ax.plot(x, y, ".", ms=5., markeredgewidth=0., color="black", alpha=.75)
                        eval_var[-1].append(data.iloc[i, :].var())
                        if abs(x) >= limit_x and y < limit_y and limit_x != None:
                            if x < 0:
                                ax.plot(x, y, ".", ms=markersize, markeredgewidth=0., color=colors[1], alpha=.4)
                            else:
                                ax.plot(x, y, ".", ms=markersize, markeredgewidth=0., color=colors[2], alpha=.4)
                        elif abs(x) >= limit_x and y > limit_y and limit_x != None:
                            if plot_ci == True:
                                ax.plot(ci,[y,y],color="black",alpha=.5, lw=.75)
                            if x < 0:
                                ax.plot(x, y, ".", ms=markersize_sign, markeredgewidth=0., color=colors[1], alpha=.9)
                                left += 1
                            else:
                                ax.plot(x, y, ".", ms=markersize_sign, markeredgewidth=0., color=colors[2], alpha=.9)
                                right += 1
                            if annotations == True and str(p) != "nan":
                                labels.append(plt.text(x, y, string.join(p.split(separator)[:max_names]), fontsize=fontsize_annotations))
                                #labels.append(" ".join(p.split(separator)[:max_names]))
                                labels_xs.append(x)
                                labels_ys.append(y)
                            #ax.text(x, y + limit_y * 0.02 + limit_y * 0.005 * random.random(),
                            #        string.join(p.split(" ")[:2], " "), horizontalalignment='center', fontsize=5)
                            above_sign += 1
                            eval_change[-1].append(abs(x))
                            eval_qvals[-1].append(y)
                        elif y > limit_y and limit_x == None:
                            ax.plot(x, y, ".", ms=markersize, markeredgewidth=0., color=colors[1], alpha=.75)
                        elif y > limit_y:
                            if x < 0:
                                ax.plot(x, y, ".", ms=markersize, markeredgewidth=0., color=colors[1], alpha=.5)
                            else:
                                ax.plot(x, y, ".", ms=markersize, markeredgewidth=0., color=colors[2], alpha=.5)
                        else:
                            ax.plot(x,y,".",ms=markersize,markeredgewidth=0.,color = colors[0], alpha=.4)
                        if y > limit_y and limit_x == None:
                            above_sign +=1
                            eval_change[-1].append(abs(x))
                            eval_qvals[-1].append(y)
                        if y > limit_y and limit_x == None and "norm" in p:
                            norm += 1
                    else:
                        if code in p:
                            ax.plot(x, y, ".", ms=markersize, markeredgewidth=0., color=colors[1], alpha=1.)
                            above_sign += 1
                        else:
                            ax.plot(x, y, ".", ms=markersize, markeredgewidth=0., color=colors[0], alpha=.1)
                    res = res.append(pd.DataFrame([[p, x, 10.**(-y), n1, n2]], columns = cols))

            #above_sign = np.round(above_sign / float(num_proteins) * 100, 1)
            eval_num.append(above_sign)
            print "... Custom above significance:", above_sign, norm
            res.to_excel(data_folder + "/SPECIFICITY/RES_expressionDE_PAULA_%s-vs-%s.xlsx"%(group1, group2), index=False)
            ax.plot([min(res.loc[:,x_col]), max(res.loc[:,x_col])],[limit_y, limit_y],"--", color = "black", lw=.5)
            if limit_x != None:
                ax.plot([limit_x, limit_x],[min(-np.log10(res.loc[:,y_col])),y_max],"--", color = "black", lw=.5)
                ax.plot([-limit_x, -limit_x],[min(-np.log10(res.loc[:,y_col])),y_max],"--", color = "black", lw=.5)
            ax.grid(False)
            if 0 < subplot < 3:
                ax.set_xlabel("fold change [glog]")
            else:
                ax.set_xlabel("fold change [log2]")
            if code == None:
                ax.text(0.05,0.95,"%s"%left, transform=ax.transAxes,ha="left", va="top", color=colors[1], alpha=.9, weight="bold")
                ax.text(0.95,0.95,"%s"%right, transform=ax.transAxes,ha="right", va="top", color=colors[2], alpha=.9, weight="bold")
                #ax.plot([], [], ".", color=colors[1], label="n=%s"%(above_sign))
                #ax.plot([], [], ".", color=colors[0], alpha=.5, label="n=%s"%(len(protein_names)-above_sign))
                #ax.legend(fancybox=False)
                #ax.set_ylim(0,)
                ax.set_xlim(min(res.loc[:,x_col]), max(res.loc[:,x_col]))
                #ax.set_ylim(min(-np.log10(res.loc[:, y_col])), max(-np.log10(res.loc[:, y_col])))

                if y_min == None:
                    ax.set_ylim(min(-np.log10(res.loc[:, y_col])), y_max)
                else:
                    ax.set_ylim(y_min, y_max)
            elif "VSN" in mode:
                ax.set_xlim(-10., 10.)
                ax.set_ylim(0, 6)
            else:
                ax.set_xlim(-1.5,1.5)
                ax.set_ylim(0, 6)
            if FDR == True:
                ax.set_ylabel("q-value [-log10]")
            else:
                ax.set_ylabel("p-value [-log10]")
    if annotations == True:
        #for label, x, y in zip(labels, labels_xs, labels_ys):
        #    #ax.text(x,y,label, ha="center", va="bottom", fontsize = fontsize_annotations, alpha=.5)
        #    ax.annotate(label,xy=(x, y), xycoords='data', xytext=(0, 10), textcoords='offset points', \
        #                ha="center", va="center", fontsize = fontsize_annotations, alpha=.7)
        adjust_text(labels, x=labels_xs, y=labels_ys, force_points=2, arrowprops=dict(arrowstyle='-',color='grey', lw=.5))
    #plt.show()
    fig.savefig(data_folder + "/SPECIFICITY/PLOT_expressionDE_PAULA_%s-vs-%s.png"%(group1, group2),dpi=300, transparent=True)