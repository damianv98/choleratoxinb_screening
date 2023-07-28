import glob, os
import statistics as st
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--proteinname', type=str, required=True)
parser.add_argument('--filenames', type=str, required=True)
parser.add_argument('--filetype', type=str, required=True)
args = parser.parse_args()
filetype = str('.')+args.filetype
protein = args.proteinname
pdblist = args.filenames.split(',')
methods = ['global', 'site', 'ddm', 'sse', 'weighted', 'sitehopper']
for m in methods:
    if not os.path.exists(m):
        os.makedirs(m)
    if m == 'sse':
        score = 'Tanimoto'
    else:
        score = 'RMSD'
    rmsdlists = []
    for i in pdblist:
        rmsdlist = []
        for j in pdblist:
            os.system('superposition -ref '+str(i)+str(filetype)+' -fit '+str(j)+str(filetype)+' -prefix '+str(m)+'/'+str(i)+str(j)+' -method '+str(m)+' -log '+str(m)+'/'+str(i)+str(j)+'.log')
            f = open(str(m)+'/'+str(i)+str(j)+'.log')
            for line in f:
                if score in line:
                    split_line = line.split(' ')
                    rmsd = split_line[3]
                    rmsd = rmsd.replace(',', '')
                    rmsdlist.append(float(rmsd))
        rmsdlists.append(rmsdlist)
    rmsdmatrix = np.stack(rmsdlists)
    df = pd.DataFrame(rmsdmatrix, index=pdblist, columns=pdblist)
    medians = []
    means = []
    for i in rmsdlists:
        medians.append(st.median(i))
        means.append(st.mean(i))
    df['median'] = medians
    df['mean'] = means
    df.index.name = m
    df.to_csv(str(protein)+' '+str(m)+'.csv', index=True, header=True, sep=' ')

writer = pd.ExcelWriter(str(protein)+' superposition results.xlsx')

for filename in glob.glob(str(protein)+'*.csv'):
    usedmethod = filename.split('\t')
    df_csv = pd.read_csv(filename, sep=' ')
    (_, f_name) = os.path.split(filename)
    (f_shortname, _) = os.path.splitext(f_name)
    df_csv.to_excel(writer, f_shortname, index=False)

writer.close()