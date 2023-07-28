import time
import random
from pathlib import Path
import re
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator

compound_df = pd.read_csv("top1000_1EEI_smiles.txt", sep=' ', header=None)
compound_df.columns = ["SMILES", "Name"]

compounds = []
for _, smiles, name in compound_df[["SMILES", "Name"]].itertuples():
    compounds.append((Chem.MolFromSmiles(smiles), name))

rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
fingerprints = [rdkit_gen.GetFingerprint(mol) for mol, idx in compounds]

print("Number of compounds converted:", len(fingerprints))
print("Fingerprint length per compound:", len(fingerprints[0]))


def tanimoto_distance_matrix(fp_list):
    """Calculate distance matrix for fingerprint list"""
    dissimilarity_matrix = []
    # Notice how we are deliberately skipping the first and last items in the list
    # because we don't need to compare them against themselves
    for i in range(1, len(fp_list)):
        # Compare the current fingerprint against all the previous ones in the list
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        # Since we need a distance matrix, calculate 1-x for every element in similarity matrix
        dissimilarity_matrix.extend([1 - x for x in similarities])
    return dissimilarity_matrix

n = len(fingerprints)
elem_triangular_matr = (n * (n - 1)) / 2
print(
 f"Elements in the triangular matrix ({elem_triangular_matr:.0f}) ==",
 f"tanimoto_distance_matrix(fingerprints) ({len(tanimoto_distance_matrix(fingerprints))})",
)

def cluster_fingerprints(fingerprints, cutoff=0.4):
    """Cluster fingerprints
    Parameters:
        fingerprints
        cutoff: threshold for the clustering
    """
    # Calculate Tanimoto distance matrix
    distance_matrix = tanimoto_distance_matrix(fingerprints)
    # Now cluster the data with the implemented Butina algorithm:
    clusters = Butina.ClusterData(distance_matrix, len(fingerprints), cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters

# Run the clustering procedure for the dataset
clusters = cluster_fingerprints(fingerprints, cutoff=0.6)
# Give a short report about the numbers of clusters and their sizes
num_clust_g1 = sum(1 for c in clusters if len(c) == 1)
num_clust_g5 = sum(1 for c in clusters if len(c) > 5)
num_clust_g25 = sum(1 for c in clusters if len(c) > 25)
num_clust_g100 = sum(1 for c in clusters if len(c) > 100)
print("total # clusters: ", len(clusters))
print("# clusters with only 1 compound: ", num_clust_g1)
print("# clusters with >5 compounds: ", num_clust_g5)
print("# clusters with >25 compounds: ", num_clust_g25)
print("# clusters with >100 compounds: ", num_clust_g100)
# NBVAL_CHECK_OUTPUT

# Plot the size of the clusters
fig, ax = plt.subplots(figsize=(15, 4))
ax.set_xlabel("Cluster index")
ax.set_ylabel("Number of molecules")
ax.bar(range(1, len(clusters) + 1), [len(c) for c in clusters], lw=5);
plt.show()
import numpy
for cutoff in numpy.arange(0.0, 1.0, 0.2):
 clusters = cluster_fingerprints(fingerprints, cutoff=cutoff)
 fig, ax = plt.subplots(figsize=(15, 4))
 ax.set_title(f"Threshold: {cutoff:3.1f}")
 ax.set_xlabel("Cluster index")
 ax.set_ylabel("Number of molecules")
 ax.bar(range(1, len(clusters) + 1), [len(c) for c in clusters], lw=5)
plt.show()

cutoff = 0.6 # Change cut-off value here
clusters = cluster_fingerprints(fingerprints, cutoff=cutoff)
# Plot the size of the clusters - save plot
fig, ax = plt.subplots(figsize=(15, 4))
ax.set_xlabel("Cluster index")
ax.set_ylabel("# molecules")
ax.bar(range(1, len(clusters) + 1), [len(c) for c in clusters])
ax.set_title(f"Threshold: {cutoff:3.1f}")
plt.show()
fig.savefig("cluster_dist_cutoff_{cutoff:4.2f}.png", dpi=300, bbox_inches="tight", transparent=True)
print(f"Number of clusters: {len(clusters)} from {len(compounds)} molecules at distance cut-off {cutoff:.2f}")
print("Number of molecules in largest cluster:", len(clusters[0]))
print(f"Similarity between two random points in same cluster: {DataStructs.TanimotoSimilarity(fingerprints[clusters[0][0]], fingerprints[clusters[0][1]]):.2f}")
print(f"Similarity between two random points in different cluster: {DataStructs.TanimotoSimilarity(fingerprints[clusters[0][0]], fingerprints[clusters[1][0]]):.2f}")

print("Twenty molecules from largest cluster:")
# Draw molecules
img = Draw.MolsToGridImage(
 [compounds[i][0] for i in clusters[1][:20]], # change the numbers for clusters[x][:y] where x is the cluster you want to see, and :y is the number of molecules you want to see from that cluster
 legends=[compounds[i][1] for i in clusters[1][:20]], # change the numbers here as well
 molsPerRow=5,
)
img.save('1EEI_cl2_20.png')


print("Twenty molecules from first 20 clusters:")
# Draw molecules
img = Draw.MolsToGridImage(
 [compounds[clusters[i][0]][0] for i in range(20)], # change range(x) where x is the number of clusters you want to see a molecule of
 legends=[compounds[clusters[i][0]][1] for i in range(20)], # change the number here as well
 molsPerRow=5,
)
img.save('output_1PZJ.png')


# Generate image
img = Draw.MolsToGridImage(
 [compounds[clusters[i][0]][0] for i in range(0, 6)],
 legends=[f"Cluster {i}" for i in range(1, 7)],
 subImgSize=(250, 250),
 useSVG=True,
)
# Patch RAW svg data: convert non-transparent to transparent background and set font size
molsvg = img.replace("opacity:1.0", "opacity:0.0").replace("12px", "20px")
# Save altered SVG data to file
with open("cluster_representatives_1EEI_6.svg", "w") as f:
 f.write(molsvg)


cluster_mols = []
for i, x in enumerate(clusters):
    cluster_mols.append([compounds[j][1] for j in clusters[i][:]])

for d, i in enumerate(cluster_mols):
    with open('cluster' + str(d+1) + '_1EEI.txt', 'w') as outfile:
        for j in i:
            outfile.write(j+"\n")