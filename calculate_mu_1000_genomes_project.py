#   calculate MU for the 1000 Genomes Project phase 3 data set
#
#   Syntax:
#   	python3 calculate_mu_1000_genomes_project.py
#
#   Output:
#       MU_result.csv
#       Each row represents a locus
#       column 1: x coordinate of the locus in the two-dimensional embedding space
#       column 2: y coordinate of the locus
#       column 3: MU value of the locus

import pandas as pd
import matplotlib.pyplot as plt
import sys # to use input parameters
import numpy as np
from sklearn import manifold

def binary_similarity(v1,v2): # v1 and v2 are discrete probability distributions
    return 1 - np.absolute(v1[0] - v2[0])

# Bhattacharyya coefficient
def Bhattacharyya_coef(v1,v2): # v1 and v2 are discrete probability distributions
    Nb = len(v1) # = 45. Number of genes 
    sim = 0.0 # similarity between v1 and v2, normalized between 0 and 1
    for i in range(Nb):
        sim += np.sqrt(v1[i]*v2[i])
    return sim

# weighted Bhattacharyya coefficient
def weighted_Bhattacharyya_coef(v1,v2): # v1 and v2 are discrete probability distributions
    Nb = len(v1) # = 45
    c = 0.5 # constant boost
    v1_weighted = np.empty(Nb)
    v2_weighted = np.empty(Nb)
    for i in range(9): # 9 columns from 0/0 to 0/8 (from 0 to 8)
        v1_weighted[i] = v1[i] * (c + i)
        v2_weighted[i] = v2[i] * (c + i)
    for i in range (9, 17): # 8 columns from 1/1 to 1/8 (from 9 to 16) 
        v1_weighted[i] = v1[i] * (c + 1 + i)
        v2_weighted[i] = v2[i] * (c + 1 + i)
    for i in range (17, 24): # 7 columns from 2/2 to 2/8 (from 17 to 23) 
        v1_weighted[i] = v1[i] * (c + 2 + i)
        v2_weighted[i] = v2[i] * (c + 2 + i)
    for i in range (24, 30): # 6 columns from 3/3 to 3/8 (from 24 to 29) 
        v1_weighted[i] = v1[i] * (c + 3 + i)
        v2_weighted[i] = v2[i] * (c + 3 + i)
    for i in range (30, 35): # 5 columns from 4/4 to 4/8 (from 30 to 34) 
        v1_weighted[i] = v1[i] * (c + 4 + i)
        v2_weighted[i] = v2[i] * (c + 4 + i)
    for i in range (35, 39): # 4 columns from 5/5 to 5/8 (from 35 to 38) 
        v1_weighted[i] = v1[i] * (c + 5 + i)
        v2_weighted[i] = v2[i] * (c + 5 + i)
    for i in range (39, 42): # 3 columns from 6/6 to 6/8 (from 39 to 41) 
        v1_weighted[i] = v1[i] * (c + 6 + i)
        v2_weighted[i] = v2[i] * (c + 6 + i)
    for i in range (42, 44): # 2 columns from 7/7 to 7/8 (from 42 to 43) 
        v1_weighted[i] = v1[i] * (c + 7 + i)
        v2_weighted[i] = v2[i] * (c + 7 + i)
    for i in range (44, 45): # 1 column from 8/8 to 8/8 (column 44) 
        v1_weighted[i] = v1[i] * (c + 8 + i)
        v2_weighted[i] = v2[i] * (c + 8 + i)
    
    normalization_factor = np.sum(v1_weighted)
    v1_weighted = v1_weighted / normalization_factor
    normalization_factor = np.sum(v2_weighted)
    v2_weighted = v2_weighted / normalization_factor
    sim = 0.0 # similarity between v1 and v2, normalized between 0 and 1
    for i in range(Nb):
        sim += np.sqrt(v1_weighted[i]*v2_weighted[i])
    return sim

# calculate distance between two similarity matrices corresponding to two genes
def Frobenius(S1, S2): 
    return np.linalg.norm(S1 - S2)

# [This section is unused.]
#
# S1 and S2 are assumed to be strictly lower diagonal
#    distance = 0.0;
#    for k in range(len(S1)): # len(S1) = Np*(N-1)/2
#        distance += (S1[k] - S2[k])**2
#    distance = np.sqrt(distance)
#    return distance


# calculate the similar matrix for a locus    
def generate_similarity_matrix(x): # x is a dataframe row, corresponding to one locus
    Np = 26 # number of populations
    Nb = 45 # = 9 * 10 / 2 = number of bins
    A = x.to_numpy()
    B = np.empty(Np*Nb)
    for i in range(Np):
        summed = np.sum(A[0, i*Nb : (i+1)*Nb])
        B[i*Nb : (i+1)*Nb] = A[0, i*Nb : (i+1)*Nb] / float(summed) # somehow we cannot reuse B

    # calculate similarity matrix
    S = np.empty(int(Np*(Np-1)/2)) # ignore the diagonal, which is always equal to 1.
    k = 0
    for i in range(Np):
        for j in range(i):
            S[k] = weighted_Bhattacharyya_coef(B[i*Nb : (i+1)*Nb], B[j*Nb : (j+1)*Nb])
            k += 1
    return S
        
# main starts here

# input data
# df = pd.read_csv('1000Genome_genotype.noX.compressed.csv', header=1, index_col = 4)
df = pd.read_csv('1000Genome_genotype.noX.compressed.withSNP.csv', header=1, index_col = 4)
# The second row in the original file is used as the header and is the population name.

NL = len(df) # Number of rows = number of loci = 57629, excluding the first two header rows.

# [Test the method only using the first NL_keep loci]
#
# NL_keep = 48000
# df.drop(df.tail(NL - NL_keep).index, inplace=True) # delete the first 7 columns and overwrite on df
# NL = NL_keep

# print(len(df.columns)) # Number of columns = 1178 including the first eight header columns.
# There are 26 populations. So, (1178 - 8)/26 = 45 = 9 * 10 / 2 columns per population
# Use the ID in column E as index. 
# Then, there are 1177 columns.

kind = df.iloc[:,3] # kind of mutation
#kind = kind.str.replace('CNV', 'OTH', case=False)
kind = kind.values.tolist() # transform to a list

for i, x in enumerate(kind):
    if x != 'DEL' and x != 'DUP' and x != 'INV' and x != 'CNV':
        kind[i]= 'OTH'

dg = df.iloc[:, range(3,5)].copy() # deep copy
df.drop(df.columns[range(7)], axis=1, inplace=True) # delete the first 7 columns and overwrite on df

Np = 26 # number of populations
SS = np.zeros((NL, int(Np*(Np-1)/2)))     
for i in range(NL):                  
# x = df.iloc[2,5:] # data row 2 (excluding the header). Index starts from 0.
    SS[i,:] = generate_similarity_matrix(df[i:i+1])

print("The similarity matrix for each locus has been created.")

if_use_Frobenius = 1
# if_use_Frobenius == 1: we use the distance matrix between loci calculated by our Frobenius distance function.
# if_use_Frobenius \neq 1: manifold.MDS calculates the distance matrix.
# The default is "if_use_Frobenius == 1"

if if_use_Frobenius == 1:
    network_distance_matrix = np.zeros((NL, NL))
    for i in range(NL):
    #    network_distance_matrix[i, i] = 0.0
        for j in range(i):
            network_distance_matrix[i, j] = Frobenius(SS[i,:], SS[j,:])
            network_distance_matrix[j, i] = network_distance_matrix[i, j]
    print("Distance matrix for locus pairs has been computed.")
# print(network_distance_matrix)
#    np.savetxt('network_distance_matrix.txt', network_distance_matrix)
#    print(np.loadtxt('sample_2.txt'))
# File input and output for numpy
# https://www.python-izm.com/data_analysis/numerical/numpy/ndarray_filerw/

    model = manifold.MDS(n_components=2, dissimilarity='precomputed', random_state=1, n_init=1)
    mds_out = model.fit_transform(network_distance_matrix)
else:
    model = manifold.MDS(n_components=2, n_init=1, max_iter=100)
    mds_out = model.fit_transform(SS)

# MDS and other clustering methods in python 
# https://qiita.com/TomHortons/items/2064870af3f3f7f2b209

x_tmp = mds_out[:,0] # string type
y_tmp = mds_out[:,1]
x = np.array(x_tmp) # somehow this is necessary. x = np.array(x) does not work
y = np.array(y_tmp)

x_ave = np.mean(x) # x_ave is approximately 0.
y_ave = np.mean(y) # y_ave is approximately 0.
MU = np.sqrt(np.square(x - x_ave) + np.square(y - y_ave)) # Euclidean distance from the origin
dg['x'] = x
dg['y'] = y
dg['MU'] = MU
dg = dg.sort_values(by='MU', ascending=False) # sort rows in the descending order of MU
dg.to_csv("MU_result.csv")
