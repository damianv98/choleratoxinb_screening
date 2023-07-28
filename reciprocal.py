hybrid_dictionary = {}
with open("D:/reciprocal_hybrid_1EEI.txt") as f:
    for line in f:
        (key, val) = line.split()
        hybrid_dictionary[key] = val

rocs_dictionary = {}
with open("D:/ROCS_1EEI_reciprocal.txt") as f:
    for line in f:
        (key, val) = line.split()
        rocs_dictionary[key] = val

hybrid_mols = list(hybrid_dictionary.keys())
rocs_mols = list(rocs_dictionary)
mols_both = set(hybrid_mols) & set(rocs_mols)

reciprocal_sum = []
for i in mols_both:
    reciprocal_sum.append(i+' '+hybrid_dictionary[i]+' '+rocs_dictionary[i]+' '+str((float(hybrid_dictionary[i])+float(rocs_dictionary[i]))))

with open('reciprocal_sum_1EEI.txt', 'w') as fp:
    fp.write('\n'.join(reciprocal_sum))