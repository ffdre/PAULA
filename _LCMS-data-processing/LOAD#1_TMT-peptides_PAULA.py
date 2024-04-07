import numpy as np
import pylab as plt
import sys, os, inspect
import pandas as pd
directory=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
file_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(file_directory)

###### settings #######
strict_quantitation = True
combine_by_mean = True
#######################

def normalize_peptides(peptides, control):
    global combine_by_mean
    peptides = peptides.div(control, axis=0)
    if combine_by_mean == True:
        peptides = peptides.mean(skipna=True)
    else:
        peptides = peptides.median(skipna=True)
    return peptides

for file in sorted(os.listdir(os.getcwd()+"/data")):
    if file.endswith("_peptides.xlsx"):
        run = file.split("_")[2]
        # print file[2:]
        print(">>> loading run #%s" % run)
        data = pd.read_excel("data/"+file)  # , index_col=concatenation_column)
        master = pd.read_excel("data/"+file[:-13]+"Master.xlsx")
        columns = data.columns
        abundance_columns = []
        sample_columns = []
        for column in columns:
            if "Abundance: " in column:
                abundance_columns.append(column)
                if "134N" not in column:
                    sample_columns.append(column)
                if "134N" in column:
                    control_column = column
        # remove non-quantified peptides
        data = data.dropna(subset=abundance_columns)
        len_1 = len(data.iloc[:,0])
        # create new column for actual number of peptides used for quantification
        master["# unique PSM"] = np.nan

        for i, row in master.iterrows():
            if np.mod(i,500) == 0:
                print("... calculating protein #%s"%i)
            protein = row["Accession"]
            if strict_quantitation == False:
                temp = data[abundance_columns].where(data["Master Protein Accessions"].str.contains(protein)).dropna(thresh=1)
            else:
                temp = data[abundance_columns].where(data["Master Protein Accessions"] == protein).dropna(thresh=1)
            temp[control_column] = temp[control_column].replace(0, np.nan)
            temp = temp.dropna(subset=[control_column]) #necessary to remove standard with intensity 0, would otherwise cause inf
            peptides = temp[sample_columns]
            control = temp[control_column]
            values = normalize_peptides(peptides, control)
            number = len(control.values)
            master.loc[master["Accession"] == protein, sample_columns] = values.values
            master.loc[master["Accession"] == protein, "# unique PSM"] = number
        master = master.drop(columns=[control_column])
        master.to_excel("data/"+file[:-13]+"Master_Normed.xlsx")


print("... all done and saved. <<<")
sys.exit()


data = data.drop_duplicates(subset=abundance_columns)
len_2 = len(data.iloc[:,0])
if len_1 != len_2:
    print("... ! duplicates detected ! # %s"%(len_1 - len_2))
print("... # of individual proteins detected: %s"%len_2)