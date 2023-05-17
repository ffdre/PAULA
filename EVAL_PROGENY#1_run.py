import pandas as pd
import numpy as np
import sys, os, inspect
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append('/Users/dressler/Documents/Core/_PythonModules')
import visualize as vs
import decoupler as dc
plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)
vs.style('science')
colors = vs.colors_indexed

try: os.mkdir("*PAULA/PROGENY")
except: pass

data = pd.read_excel("*PAULA/NORMICS_results/RES_data_normalizedNORMICSmed-log2.xlsx")
index = pd.merge(
    pd.read_excel("*PAULA/NORMICS_results/RES_data_proteins.xlsx"),
    pd.read_excel("*PAULA/*RES_data_functional.xlsx")[["Accession","Gene Symbol"]],
    how='left', on="Accession")
data = data.set_index(index["Gene Symbol"].values)
data = 2.**data
#data["Index"] = index["Gene Symbol"].values
#data = data.replace(np.nan, 0.)
#data = data.astype(float)

print(data)
print(data.index.values)

model = dc.get_progeny(organism='human', top=100)
#intersect = pd.merge(data["Index"],model, left_on="Index", right_on="target", how='inner')
#intersect.to_excel("*PAULA/PROGENY/Ctrl.xlsx")


print(model)

print(data.values)
print(list(data.columns))
print(list(data.index.values))
print(len(data.values[0]))
print(len(data.values[1]))
estimate, pvals = dc.run_mlm(mat=[data.T.values, list(data.columns),list(data.index.values)], net=model, verbose=True, min_n=3) #source='source', target='target', weight='weight', min_n=1

print(estimate)
print(pvals)
#data.to_excel("*PAULA/PROGENY/*RES_progeny.xlsx")
estimate.to_excel("*PAULA/PROGENY/*RES_progeny_estimate.xlsx")
pvals.to_excel("*PAULA/PROGENY/*RES_progeny_pvals.xlsx")







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
