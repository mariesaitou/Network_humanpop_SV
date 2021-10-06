#   calculate MU for the simulated data, by specifiying the input and output files as arguments 
#   Syntax
#
#   python3 calculate_mu_simulated_input_output_specified.py
#
#  input file: SLiM.input.1.csv in_file out_file
# 
# column A: locus (3000 in total)
# column B-AK: African population (36 genotypes, 100 individuals)
# column AL-BU: European population (36 genotypes, 100 individuals)
# column BV-DE: Asian population (36 genotypes, 100 individuals)
# column DF: mutation type 
#
# output format
#
# locus, x, y, MU, and "type" in the csv file. 

import pandas as pd
import matplotlib.pyplot as plt
import sys # to use input arguments
import numpy as np
from sklearn import manifold

def binary_similarity(v1,v2): # v1 and v2 are discrete probability distributions
    return 1 - np.absolute(v1[0] - v2[0])

def Bhattacharyya_coef(v1,v2): # v1 and v2 are discrete probability distributions
    Nb = len(v1) # = 36
    sim = 0.0 # similarity between v1 and v2, normalized between 0 and 1
    for i in range(Nb):
        sim += np.sqrt(v1[i]*v2[i])
    return sim

def weighted_Bhattacharyya_coef(v1,v2): # v1 and v2 are discrete probability distributions
    Nb = len(v1) # = 36 genotypes
    c = 0.5 # constant boost
    v1_weighted = np.empty(Nb)
    v2_weighted = np.empty(Nb)
    for i in range(8): # 8 columns from 1/1 to 1/8 (from 0 to 7)
        v1_weighted[i] = v1[i] * (c + i)
        v2_weighted[i] = v2[i] * (c + i)
    for i in range (8, 15): # 7 columns from 2/2 to 2/8 (from 8 to 14) 
        v1_weighted[i] = v1[i] * (c + 1 + i)
        v2_weighted[i] = v2[i] * (c + 1 + i)
    for i in range (15, 21): # 6 columns from 3/3 to 3/8 (from 15 to 20) 
        v1_weighted[i] = v1[i] * (c + 2 + i)
        v2_weighted[i] = v2[i] * (c + 2 + i)
    for i in range (21, 26): # 5 columns from 4/4 to 4/8 (from 21 to 25) 
        v1_weighted[i] = v1[i] * (c + 3 + i)
        v2_weighted[i] = v2[i] * (c + 3 + i)
    for i in range (26, 30): # 4 columns from 5/5 to 5/8 (from 26 to 29) 
        v1_weighted[i] = v1[i] * (c + 4 + i)
        v2_weighted[i] = v2[i] * (c + 4 + i)
    for i in range (30, 33): # 3 columns from 6/6 to 6/8 (from 30 to 32) 
        v1_weighted[i] = v1[i] * (c + 5 + i)
        v2_weighted[i] = v2[i] * (c + 5 + i)
    for i in range (33, 35): # 2 columns from 7/7 to 7/8 (from 33 to 34) 
        v1_weighted[i] = v1[i] * (c + 6 + i)
        v2_weighted[i] = v2[i] * (c + 6 + i)
    for i in range (35, 36): # 1 column 8/8 (column 35) 
        v1_weighted[i] = v1[i] * (c + 7 + i)
        v2_weighted[i] = v2[i] * (c + 7 + i)
    
    normalization_factor = np.sum(v1_weighted)
    v1_weighted = v1_weighted / normalization_factor
    normalization_factor = np.sum(v2_weighted)
    v2_weighted = v2_weighted / normalization_factor
    sim = 0.0 # similarity between v1 and v2, normalized between 0 and 1
    for i in range(Nb):
        sim += np.sqrt(v1_weighted[i]*v2_weighted[i])
    return sim

def Frobenius(S1, S2): # calculate distance between two similarity matrices
    return np.linalg.norm(S1 - S2)

# S1 and S2 are assumed to be strictly lower diagonal
#    distance = 0.0;
#    for k in range(len(S1)): # len(S1) = Np*(N-1)/2
#        distance += (S1[k] - S2[k])**2
#    distance = np.sqrt(distance)
#    return distance

#    return np.sqrt(np.dot(S1, S1) - 2*np.dot(S1, S2) + np.dot(S2, S2))
   
#    return np.sqrt(np.sum((S1 - S2)**2)) 
    
def generate_similarity_matrix(x): # x is a dataframe row, corresponding to one locus
    Np = 3 # number of populations
    Nb = 36 # = 8 * 9 / 2 = number of bins
#    print(x.columns[5:(len(x.columns)):Nb])
#    A = x.iloc[:,7:(len(x.columns))].to_numpy()
    A = x.to_numpy()
#     print(A.shape)
    B = np.empty(Np*Nb)
#    print(B.shape)
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
        

argc = len(sys.argv)
if argc != 3:
    sys.exit('Usage: python calculated_mu_simulated_input_output_specified.py input_file_name output_file_name')
    
    
# columns 1 and 2: two nodes to form an event
# column 3: timestamp of the event

# df = pd.read_csv('1000Genome_genotype.noX.compressed.csv', header=1, index_col = 4)

# df = pd.read_csv('SLiM.input.1.csv', header=0, index_col = 0)
df = pd.read_csv(sys.argv[1], header=0, index_col = 0)

NL = len(df) # rows = 57629 excluding the first two header rows.
# The second row in the original file is used as the header and is the population name.

#NL_keep = 48000
#df.drop(df.tail(NL - NL_keep).index, inplace=True) # delete the first 7 columns and overwrite on df
#NL = NL_keep

# print(len(df.columns)) # Number of columns = 109 excluding the first column of the data file (= locus) but including the last column (= mutation type).
# There are 3 populations and 36 genotypes. So, 3*36 = 108 data columns.
# The first column of the data = locus.
# The last (i.e., 110th) column of the data = mutation type.

mutation_type = df.iloc[:,108] # mutation type
# print(type(mutation_type)) # type = panda series
# mutation_type = mutation_type.str.replace('CNV', 'OTH', case=False)
mutation_type = mutation_type.values.tolist() # transform to a list

# print(mutation_type)

dg = df.iloc[:, range(108, 109)].copy() # deep copy
df.drop(df.columns[108], axis=1, inplace=True) # delete the last column and overwrite on df
# print(len(df.columns))

Np = 3 # number of populations
SS = np.zeros((NL, int(Np*(Np-1)/2)))     
for i in range(NL):                  
# x = df.iloc[2,5:] # data row 2 (excluding the header). Index starts from 0.
    SS[i,:] = generate_similarity_matrix(df[i:i+1])

print("The similarity matrix for each locus has been created.")

if_use_Frobenius = 1

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


# labels = np.genfromtxt("park_time.csv",delimiter=",",usecols=0,dtype=str)

#for label, x, y in zip(labels, pos[:, 0], pos[:, 1]):
#    plt.annotate(
#        label,
#        xy = (x, y), xytext = (70, -20),
#        textcoords = 'offset points', ha = 'right', va = 'bottom',
#        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
#        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')
#    )
#
## https://dev.classmethod.jp/statistics/mds-for-parks/

x_tmp = mds_out[:,0] # string type
y_tmp = mds_out[:,1]
x = np.array(x_tmp) # somehow this is necessary. x = np.array(x) does not work
y = np.array(y_tmp)

x_ave = np.mean(x)
y_ave = np.mean(y)
# print("center = (", x_ave, y_ave, ")") # x_ave, y_ave \approx 0
d_from_ave = np.sqrt(np.square(x - x_ave) + np.square(y - y_ave))
#df = df.drop(range(5, len(df.columns)), axis=1)
dg['x'] = x
dg['y'] = y
dg['d_from_ave'] = d_from_ave
dg = dg.sort_values(by='d_from_ave', ascending=False)
print(dg)
dg = dg[['x','y','d_from_ave','type']]
print(dg)
# print(df)
# dg.to_csv("MDS_simulated_result.csv")
dg.to_csv(sys.argv[2])

# pd.options.mode.chained_assignment = None  # default='warn'