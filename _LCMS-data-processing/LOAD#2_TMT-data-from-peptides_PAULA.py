### import modules
import urllib
import numpy as np
import pylab as plt
import pandas as pd
import sys, os, inspect
import scipy
import string
from operator import itemgetter

plt.close('all')
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

################################
normalize_batch = False #has been normalized by peptide normalization, this function is thus deprecated
################################

instructions = pd.read_excel("*columns.xlsx")
columns_to_remove = list(instructions["columns_to_remove"].dropna())
columns_to_rename = list(instructions["columns_to_rename"].dropna())
columns_to_keep = list(instructions["columns_to_keep"].dropna())
#concatenation_column = list(instructions["column_for_concat"])[0]

samples = pd.read_excel("**CTRL_220111_Luebeck_TMT-sets_curated_STRICT-FINAL.xlsx")
res = pd.DataFrame()

i = 0
for file in sorted(os.listdir(os.getcwd()+"/data")):
    if file.endswith("_Master_Normed.xlsx") and "*" not in file:
        run = file.split("_")[2]
        #print file[2:]
        print(">>> loading run #%s"%run)
        try:
            data = pd.read_excel("data/"+file).iloc[:,1:]#, index_col=concatenation_column)
            columns = data.columns

            for column in columns:
                if column in columns_to_rename:
                    data = data.rename(columns={column: run + "_" + column})
                else:
                    for column_to_remove in columns_to_remove:
                        if column_to_remove in column:
                            data = data.drop(columns = [column])
                if "Abundance" in column:
                    if "134N" not in column:
                        sample = samples.loc[samples["LCMS"] == column, "ID_set"].values[0]
                        #print sample
                        data = data.rename(columns = {column: sample})
                    if "134N" in column:
                        data = data.drop(columns=[column])


            #data.to_excel("temp.xlsx")
            print("... concatenating run #%s" % run)
            if i == 0:
                res = data
            else:
                res = pd.merge(res, data, how='outer', on=columns_to_keep)# left_index=True, right_index=True)#axis=1)
            i += 1

            del data
            del columns

            ### reorder columns
            A, B, C = [], [], []
            for c in res.columns:
                if "U-" in c:
                    C.append(c)
                elif c[0] == "H" and c[3]=="_":
                    B.append(c)
                elif c[0] == "L" and c[3]=="_":
                    B.append(c)
                else:
                    A.append(c)
            res = res[A+sorted(B, key=lambda x: x[4:]) + C]  #itemgetter(4))+C]
        except:
            print("... failed to load %s"%file)
            pass


res.to_excel("*RES_data_raw_normalize_batch=TRUE.xlsx", index=False)
#res.to_csv("*RES_data_raw_normalize_batch=TRUE.csv", index=False)

print("... all done and saved. <<<")