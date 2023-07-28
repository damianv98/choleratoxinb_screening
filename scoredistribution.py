import seaborn as sns
import pandas as pd
import matplotlib as plt
moldata = pd.read_csv('D:/ROCS_1EEI_nodups.txt', sep=" ", header=None, names=["Molname", "TanimotoCombo Score"])
actives = pd.read_csv('D:/1EEI_actives_rocs.txt', sep=" ", header=None, names=["Molname", "TanimotoCombo Score"])
displot = sns.histplot(moldata["TanimotoCombo Score"], bins=200)
for i in actives["TanimotoCombo Score"]:
     displot.axvline(x=i, color='red', alpha=0.4)
displot.set(ylabel="Number of molecules")
displot.legend(labels=["Actives"])