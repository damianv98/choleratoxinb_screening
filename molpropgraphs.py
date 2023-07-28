import matplotlib.pyplot as plt
import seaborn as sbn
import dask.dataframe as dd
import pandas as pd
from pandas.api.types import is_numeric_dtype

actives = dd.read_csv("actives.txt", header=0, sep='\t')
database = dd.read_csv('database.txt', header=0, sep='\t')
properties = []
for col in actives.columns:
    properties.append(col)
properties.remove('SMILES')
properties.remove('Name')
def propertyplots(properties):
    for i in properties[:20]:
        plot = sbn.displot(database[i], bins=50)
        plot.set(xlabel=i)
        for j in actives[i]:
            plt.axvline(j, color='red')
        if (is_numeric_dtype(database[i].all)):
            plot.set(xlim=min((database[i]), max(database[i])))
        filename = i+'.png'
        plot.savefig(filename)
pd.Categorical
propertyplots(properties)