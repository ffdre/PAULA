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
#from lifelines.utils import restricted_mean_survival_time
from lifelines.utils import median_survival_times,restricted_mean_survival_time
import string
import matplotlib as mpl

vs.style('science') #, n=10
colors = vs.colors_indexed

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)




### set basic parameters ###
data_folder = "*PAULA"
k_chosen = input("k = ")
k_chosen = int(k_chosen)
fig = plt.figure("PLOT", [16,12])
fig.subplots_adjust(hspace=.4, wspace=.6)
for c,p,l in zip(["<3","==1","==2",".isin([3,4,5])","==3","==4","==5",">5","==6","==7","==8"],[1,6,11,2,7,12,17,3,8,13,18],[0,0,1,0,0,0,1,0,0,0,1]):
    condition = c
    combine = "n"
    position = p
    xlab = l
    max_follow = 60
    legend_loc = 4
    cluster_names = ["PAULA I", "PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III"]
    cmap = "Spectral"
    restricted = ""#_restricted"
    colors = vs.colors_cluster
    ############################

    data_protein_names = pd.read_excel(data_folder + "/CLUSTER/RES_data_proteins.xlsx")
    num_proteins = len(data_protein_names.iloc[:, 0])
    data_clinical = pd.read_excel(data_folder + "/*PAULA_clinical_NMFcluster.xlsx")
    descriptives = pd.read_excel(data_folder + "/CLUSTER/figures/descriptives.xlsx")
    data_survival = pd.read_excel(data_folder + "/*PAULA_survival.xlsx")

    if combine == "y":
        data_clinical["cluster k=%s"%k_chosen] = data_clinical["cluster k=%s"%k_chosen].replace([1,2,3,4],[1,1,1,2])

    clusters = np.sort(data_clinical["cluster k=%s"%k_chosen].unique())
    #### subset ####
    exec("data_clinical = data_clinical.where(data_clinical['Stage'] %s).dropna(thresh=1)"%condition)
    ################

    data_clinical = pd.merge(data_clinical, data_survival, how="inner", on="ID_set")
    N = len(data_clinical.iloc[:,0])


    #gs = mpl.gridspec.GridSpec(5,5,hspace=.4, wspace=.3)

    #ax = fig.add_subplot(gs[position,-1])
    ax = fig.add_subplot(5,5,position)

    Sum = 0
    res = pd.DataFrame()
    res["cluster"] = np.array([np.nan]*k_chosen)
    res["n="] = np.array([np.nan]*k_chosen)
    res["%"] = np.array([np.nan]*k_chosen)
    res["median survival"] = np.array([np.nan]*k_chosen)

    res["ci_low"] = np.array([np.nan]*k_chosen)
    res["ci_high"] = np.array([np.nan]*k_chosen)
    res["res. mean survival"] = np.array([np.nan]*k_chosen)
    #res["res. mean survival_low"] = np.array([np.nan]*k_chosen)
    #res["res. mean survival_high"] = np.array([np.nan]*k_chosen)
    for i,cluster in enumerate(clusters):
        temp = data_clinical.where(data_clinical["cluster k=%s"%k_chosen] == cluster).dropna(thresh=1)
        durations = temp["OS%s"%restricted]
        event_observed = temp["Event%s"%restricted]
        number = len(durations)

        km = KaplanMeierFitter()
        res.at[i, "cluster"] = cluster
        res.at[i, "n="] = number
        res.at[i, "%"] = number/float(N)*100.
        try:


            Fit = km.fit(durations, event_observed, label="_nolegend_")
            res.at[i, "median survival"] = Fit.median_survival_time_
            res.at[i, "res. mean survival"] = restricted_mean_survival_time(Fit,t=24.)
            res.at[i, "ci_low"] = median_survival_times(Fit.confidence_interval_).iloc[0,0]
            res.at[i, "ci_high"] = median_survival_times(Fit.confidence_interval_).iloc[0,1]
            #res.at[i, "ci_high"] = median_survival_times(Fit.confidence_interval_)[1]

            #print(Fit.survival_function_)
            print("...median survival:", Fit.median_survival_time_)
            print("ci", median_survival_times(Fit.confidence_interval_))
            #print("...restricted_mean_survival:",restricted_mean_survival_time(Fit, t=60))
            #Fit.event_table.to_excel(data_folder + "/CLUSTER/figures/RES_survival-_k%s_cond%s_cluster=%s_event-table.xlsx"%(k_chosen,condition,cluster), index=False)

            Fit.plot(ax = ax, ci_show = False, ci_alpha=.05, show_censors = True, color = colors[i])
            if "V" in cluster_names[i]:
                filler = " "*(14-len(cluster_names[i]))
            else:
                filler = " "*(15-len(cluster_names[i]))
            #ax.plot([None], [0.], ".", marker="s", markeredgewidth=.0, color=colors[i], label = '%s %s(n=%s)' % (cluster_names[i],"", number))
            #ax.plot([None], [0.], ".", marker="s", markeredgewidth=.0, color=colors[i], label = 'n=%s' % (number))
            ax.plot([max_follow-0.95*max_follow], [0.45-0.07*(i)], ".", marker="s", markeredgewidth=.0, color=colors[i])
            ax.text(max_follow-0.9*max_follow, 0.45-0.07*(i), 'n=%s' % (number), fontsize=7, va='center', ha='left')
        except:
            print("... skipped cluster #%s with #%s patients"%(i,len(durations)))
            pass
    #res.to_excel(data_folder + "/CLUSTER/figures/RES_survival-_k%s_cond%s.xlsx"%(k_chosen,condition), index=False)
    print(colors)
    ax.get_legend().remove()
    #ax.legend(fancybox=False, loc=legend_loc)
    ax.grid(False)
    ax.axhline(y=.5, linestyle="--", lw=1., color="black", alpha=.5)
    #ax.plot([0, max_follow], [.5, .5], "--", color = "grey", lw=1.)

    #ax.legend(fancybox=False)
    print("N=%s"%Sum)
    import statistics as ss


    print("... logrank")
    p_value = lstat.multivariate_logrank_test(data_clinical["OS_restricted"], data_clinical["Event_restricted"], data_clinical["cluster k=%s"%k_chosen]).summary
    p_value = p_value.loc[0,"p"]
    from decimal import Decimal

    if p_value < 0.0001:
        ax.text(0.05, 0.05, "p=%.2E" % Decimal(p_value), transform=ax.transAxes, size=7, weight='regular', ha='left', va='bottom')
    else:
        ax.text(0.05, 0.05, "p=%s" % np.round(p_value,4), transform=ax.transAxes, size=7, weight='regular', ha='left', va='bottom')
    ax.set_xlim(0,max_follow)
    ax.set_ylim(0,1.05)

    if xlab==1:
        ax.set_xlabel("Overall survival [months]")
    else:
        ax.set_xlabel("")
    ax.set_ylabel("Proportion surviving [1]")

fig.savefig(data_folder + "/CLUSTER/figures/PLOT_PUBL_survival_k%s_cond=ALL.png"%(k_chosen),dpi=300, transparent=True)
plt.show()

for i,c in enumerate(np.sort(data_clinical["cluster k=%s"%k_chosen].dropna().unique())):
    for j in range(i+1,len(data_clinical["cluster k=%s"%k_chosen].unique())):
        print("")
        print("... cluster #%s vs. cluster #%s" %(i+1,j+1))
        #try:
            #ss.deltaquant(os[i],os[j],names = ["cluster #%s" %(i),"cluster #%s" %(j)])
        print("... logrank")
        lstat.logrank_test(data_clinical["OS_restricted"].where(data_clinical["cluster k=%s"%k_chosen] == i).dropna(), \
                           data_clinical["OS_restricted"].where(data_clinical["cluster k=%s" % k_chosen] == j).dropna(), \
                           data_clinical["Event_restricted"].where(data_clinical["cluster k=%s"%k_chosen] == i).dropna(), \
                           data_clinical["Event_restricted"].where(data_clinical["cluster k=%s" % k_chosen] == j).dropna()).print_summary()
        #except: pass