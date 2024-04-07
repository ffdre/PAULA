import pandas as pd
import numpy as np
import sys, os, inspect
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
sys.path.append(file_directory + '/_modules')
import visualize as vs
import matplotlib as mpl
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
vs.style('science')
colors = vs.colors_indexed

data_folder = "*PAULA/SPECIFICITY"

fig = plt.figure("PLOT", [12,12])
#fig.subplots_adjust(hspace=.4, wspace=.3)
gs = mpl.gridspec.GridSpec(10,4,hspace=.4,wspace=.3, width_ratios=[1.,1.23,0.77,1.])

UP, UP_rel = [], []
DN, DN_rel = [], []
UP_fc, UP_fc_rel = [], []
DN_fc, DN_fc_rel = [], []
for file in sorted(os.listdir(os.getcwd()+"/%s"%data_folder)):
    if file.endswith(".xlsx") and "*" not in file and "pTNM" in file:
        file_id = file.split("_")[-1]
        print(">>> loading file #%s"%file_id)
        #cols = ["protein", x_col, y_col, "N1", "N2"]
        data = pd.read_excel(data_folder + "/%s"%file)
        data = data.where(data.iloc[:,2] <= .05).dropna()
        UP.append( len(data.iloc[:,1].where(data.iloc[:,1] > 0).dropna()) )
        UP_rel.append(len(data.iloc[:, 1].where(data.iloc[:, 1] > 0.5).dropna()))
        DN.append( len(data.iloc[:,1].where(data.iloc[:,1] < 0).dropna()) )
        DN_rel.append(len(data.iloc[:, 1].where(data.iloc[:, 1] < -0.5).dropna()))
        UP_fc.append( data.iloc[:,1].where(data.iloc[:,1] > 0).dropna().mean() )
        UP_fc_rel.append( data.iloc[:,1].where(data.iloc[:,1] > 0.5).dropna().mean() )
        DN_fc.append( data.iloc[:,1].where(data.iloc[:,1] < 0).dropna().mean() )
        DN_fc_rel.append( data.iloc[:,1].where(data.iloc[:,1] < -0.5).dropna().mean() )

print UP
groups_stage = ["pTa low", "pTa high", "pTis", "pT1a/b/G1/G2", "pT1c/G3", "pT2", ">pT2","pN1"]
groups_cluster = ["PAULA I", "PAULA IIa", "PAULA IIb", "PAULA IIc", "PAULA III"]

ax = fig.add_subplot(gs[0,1])
ax.set_title("pTNM", fontsize=10.)
x = range(len(UP))
ax.bar(x,UP,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,UP_rel,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,-1*np.array(DN),color=colors[0],alpha=.75, lw=.0)
ax.bar(x,-1*np.array(DN_rel),color=colors[0],alpha=.75,  lw=.0)
ax.axhline(y=0,color="black", lw=1.)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', top=False, bottom=False, labelbottom=False)
ax.tick_params(axis='y', right=False)
ax.set_yticks([500,0.,-500])
ax.set_yticklabels([500,0.,500])
ax.set_xticks(x)
ax.set_xlim(-1,len(x))
ax.set_xticklabels(groups_stage,rotation=45.,ha='right', fontsize=8.)
ax.grid(False)
ax.set_ylabel("Sig. DE [1]")

ax = fig.add_subplot(gs[1,1])
x = range(len(UP))
ax.bar(x,UP_fc,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,UP_fc_rel,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,np.array(DN_fc),color=colors[0],alpha=.75, lw=.0)
ax.bar(x,np.array(DN_fc_rel),color=colors[0],alpha=.75,  lw=.0)
ax.axhline(y=0,color="black", lw=1.)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', top=False, bottom=False)
ax.tick_params(axis='y', right=False)
ax.set_yticks([.5,0.,-.5])
ax.set_yticklabels([.5,0.,.5])
ax.set_xticks(x)
ax.set_xlim(-1,len(x))
ax.set_xticklabels(groups_stage,rotation=45.,ha='right', fontsize=8.)
ax.grid(False)
ax.set_ylabel("Median FC [log$_2$]")



UP, UP_rel = [], []
DN, DN_rel = [], []
UP_fc, UP_fc_rel = [], []
DN_fc, DN_fc_rel = [], []
for file in sorted(os.listdir(os.getcwd()+"/%s"%data_folder)):
    if file.endswith(".xlsx") and "*" not in file and "PAULA" in file:
        file_id = file.split("_")[-1]
        print(">>> loading file #%s"%file_id)
        #cols = ["protein", x_col, y_col, "N1", "N2"]
        data = pd.read_excel(data_folder + "/%s"%file)
        data = data.where(data.iloc[:,2] <= .05).dropna()
        UP.append( len(data.iloc[:,1].where(data.iloc[:,1] > 0).dropna()) )
        UP_rel.append(len(data.iloc[:, 1].where(data.iloc[:, 1] > 0.5).dropna()))
        DN.append( len(data.iloc[:,1].where(data.iloc[:,1] < 0).dropna()) )
        DN_rel.append(len(data.iloc[:, 1].where(data.iloc[:, 1] < -0.5).dropna()))
        UP_fc.append( data.iloc[:,1].where(data.iloc[:,1] > 0).dropna().mean() )
        UP_fc_rel.append( data.iloc[:,1].where(data.iloc[:,1] > 0.5).dropna().mean() )
        DN_fc.append( data.iloc[:,1].where(data.iloc[:,1] < 0).dropna().mean() )
        DN_fc_rel.append( data.iloc[:,1].where(data.iloc[:,1] < -0.5).dropna().mean() )

ax = fig.add_subplot(gs[0,2])
ax.set_title("PAULA", fontsize=10.)
x = range(len(UP))
ax.bar(x,UP,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,UP_rel,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,-1*np.array(DN),color=colors[0],alpha=.75, lw=.0)
ax.bar(x,-1*np.array(DN_rel),color=colors[0],alpha=.75,  lw=.0)
ax.axhline(y=0,color="black", lw=1.)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', top=False, bottom=False, labelbottom=False)
ax.tick_params(axis='y', right=False)
ax.set_yticks([500,0.,-500])
ax.set_yticklabels([500,0.,500])
ax.set_xticks(x)
ax.set_xlim(-1,len(x))
ax.set_xticklabels(groups_cluster,rotation=45.,ha='right', fontsize=8.)
ax.grid(False)
ax.set_ylabel("Sig. DE [1]")

ax = fig.add_subplot(gs[1,2])
x = range(len(UP))
ax.bar(x,UP_fc,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,UP_fc_rel,color=colors[1],alpha=.75, lw=.0)
ax.bar(x,np.array(DN_fc),color=colors[0],alpha=.75, lw=.0)
ax.bar(x,np.array(DN_fc_rel),color=colors[0],alpha=.75,  lw=.0)
ax.axhline(y=0,color="black", lw=1.)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', top=False, bottom=False)
ax.tick_params(axis='y', right=False)
ax.set_yticks([.5,0.,-.5])
ax.set_yticklabels([.5,0.,.5])
ax.set_xticks(x)
ax.set_xlim(-1,len(x))
ax.set_xticklabels(groups_cluster,rotation=45.,ha='right', fontsize=8.)
ax.grid(False)
ax.set_ylabel("Median FC [log$_2$]")

fig.savefig(data_folder +"/PLOT_DE-across-groups.png", dpi=300, transparent = True)
plt.show()
sys.exit()



estimate = pd.read_excel("*PAULA/PROGENY/*RES_progeny_estimate.xlsx", index_col=0)
pathways = estimate.columns[1:]
print pathways, estimate.columns
pvals = pd.read_excel("*PAULA/PROGENY/*RES_progeny_pvals.xlsx", index_col=0)
data_clinical = pd.read_excel("*PAULA/*PAULA_clinical_NMFcluster.xlsx")

if 0 == 1:
    fig0 = plt.figure("CTRL")
    ax = fig0.add_subplot(111)
    for i, D in estimate.iterrows():
        for DD,j  in enumerate(estimate.columns):
            #print i,j
            ax.plot([estimate.loc[i,j]],[pvals.loc[i,j]],".", color="black")
    plt.show()

estimate = pd.merge(estimate, data_clinical, how='left', left_index=True, right_on="ID_set")
pvals = pd.merge(pvals, data_clinical, how='left', left_index=True, right_on="ID_set")



fig = plt.figure("PLOT", [12,12])
gs = mpl.gridspec.GridSpec(12, 6, hspace=.4, wspace=.1, width_ratios = [2.,9.,5.,.25,5.,20])
#vmax = np.max([np.abs(estimate[pathways].min().min()), estimate[pathways].max().max()])

## subplot1
ax = fig.add_subplot(gs[-2:, 0])
ax.set_title("All", fontsize=10.)
data = estimate.groupby("binary_group_number")[pathways].median()
#data = np.log(data.div(data.loc[0,:],axis="columns"))
vmax = data.max().max()#np.max([np.abs(data.min().min()), data.max().max()])
vmax = 2.
ax.imshow(data.T,cmap="RdBu_r",aspect="auto", vmin=-vmax, vmax=vmax)
ax.grid(False)
ax.set_yticks(range(len(pathways)))
ax.set_yticklabels(pathways, fontsize=8.)
ax.set_xticks(range(len(data.index)))
ax.set_xticklabels(["Healthy", "Tumor"], rotation =45., ha='right', fontsize=8.)
ax.tick_params(axis='x',bottom=False,top=False)
ax.tick_params(axis='y',left=False,right=False)

## subplot2
ax = fig.add_subplot(gs[-2:, 1])
ax.set_title("pTNM", fontsize=10.)
data = estimate.groupby("Stage")[pathways].median()
#data = np.log(data.div(data.loc[0,:],axis="columns"))
#vmax = data.max().max()#np.max([np.abs(data.min().min()), data.max().max()])
print vmax
ax.imshow(data.T,cmap="RdBu_r",aspect="auto", vmin=-vmax, vmax=vmax)
ax.grid(False)
ax.set_yticks(range(len(pathways)))
ax.set_yticklabels(pathways)
ax.set_xticks(range(len(data.index)))
ax.set_xticklabels(groups_stage, rotation =45., ha='right', fontsize=8.)
ax.tick_params(axis='x',bottom=False,top=False)
ax.tick_params(axis='y',left=False,right=False, labelleft=False)

## subplot2
ax = fig.add_subplot(gs[-2:, 2])
ax.set_title("PAULA", fontsize=10.)
data = estimate.groupby("cluster k=5")[pathways].median()
#data = np.log(data.div(data.loc[0,:],axis="columns"))
#vmax = data.max().max()#np.max([np.abs(data.min().min()), data.max().max()])
print vmax
ax.imshow(data.T,cmap="RdBu_r",aspect="auto", vmin=-vmax, vmax=vmax)
ax.grid(False)
ax.set_yticks(range(len(pathways)))
ax.set_yticklabels(pathways)
ax.set_xticks(range(len(data.index)))
ax.set_xticklabels(groups_cluster, rotation =45., ha='right', fontsize=8.)
ax.tick_params(axis='x',bottom=False,top=False)
ax.tick_params(axis='y',left=False,right=False, labelleft=False)

cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax = fig.add_subplot(gs[-1,3])
ax.grid(False)
#ax.axis('off')
#ax.set_title("", rotation=45., fontsize=8., ha='left', va='bottom')
cb = fig.colorbar(sm,aspect=15., fraction=1.,
             cax=ax, orientation='vertical')
cb.set_label(label='Activity\nscore [1]',weight='regular', fontsize=7.)
cb.set_ticks(np.linspace(-vmax,vmax,3))
ax.tick_params(labelsize=7.)

#plt.colorbar(plt.cm.ScalarMappable(norm=data, cmap="RdBu_r"),ax=ax, fraction=.1)

fig.savefig("*PAULA/PROGENY/PLOT_progeny.png", dpi=300, transparent=True)
plt.show()












sys.exit()
########################################### PLOT ###########################################
data = pd.read_excel("data.xlsx")
fig = plt.figure("PLOT",[8,8])
ax = fig.add_subplot(221)
ax.bar(0,data.iloc[0,1], color = colors[0], alpha=.5)
for i in range(1,6):
    ax.bar(i,data.iloc[i,1], color = colors[1], alpha=.5)
ax.set_xticks(range(6))
ax.set_xticklabels(data.iloc[:,0], rotation=45., ha='left', va='bottom')
ax.set_ylabel("Number of proteins identified [1]")
ax.tick_params(axis="x", bottom=True, top=True, labelbottom=False, labeltop=True)
ax.grid(False)
fig.savefig("PLOT_data_reconstitution#49.png", dpi=300)
plt.show()


sys.exit()
########################################### Pandas Quick Reference Manual ###########################################
### load data
data = pd.read_excel("some_file.xlsx")
data = pd.read_excel("some_file.xlsx", index_col = 0) # ! can be convenient to set index_col to identifier

### load more data
for file in sorted(os.listdir(os.getcwd()+"/data")):
    if file.endswith(".xlsx") and "*" not in file:
        file_id = file.split("_")[0]
        print(">>> loading file #%s"%file_id)

### change value in pandas dataframe:
data.at[i,"some column"] = new_value
data.iat[i,j] = new_value
data.loc[i,"some column"] = new_values
data.iloc[i,j] = new_values

### filter for values
subset = data.where(data["some column"] == some_value).dropna(thresh=1) # ! important as otherwise subset rows with missing values are lost too!
# AND
subset = data.where((data["some column"] == some_value) & (data["some column"] == some_value)).dropna(thresh=1) # ! important as otherwise subset rows with missing values are lost too!
# OR
subset = data.where((data["some column"] == some_value) | (data["some column"] == some_value)).dropna(thresh=1) # ! important as otherwise subset rows with missing values are lost too!
# list
subset = data[data['some column'].isin([3, 6])]

### find specific value
specific_value = data["specific column"]

### find specific row
specific_row = data.loc[data.loc[:,"some column"] == some_value].iloc[0, :] # first hit

### get only array not pd.Series
array = data.values

### get unique values, sorted
uniques = np.sort(data["some_column"].unique())

### get boolean array of columns
sample_columns = ~data.columns.isin(['sort', 'accession'])

### sort data
data = data.sort_values(by="sort", ascending=False)
data = data.reset_index(drop=True) # ! important when iterrated

### merge data
new_data = pd.merge(data_left, data_right, how="inner/left/right/outer", on="some_column or list of columns", \
                    left_on="id column left", right_on="id column right")

### save data
data.to_excel("file_name.xlsx", index=False) # ! avoid unnecessary saving, consumes a lot of memory and slows script down

### ensure freeing memory
del data
del all_other_references_to_data
