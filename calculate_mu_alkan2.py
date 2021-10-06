# Calculate MU from the similarity matrix for the data provided by Can Alkan
# 
# Take similarity_matrix_pop_to_pop.csv as input
#   Syntax
#
#   python3 calculate_mu_alkan.py
#

import pandas as pd
import gc # garbage collection (to release memory)
import matplotlib.pyplot as plt
import sys # to use input parameters
import numpy as np
from sklearn import manifold

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
    
param = sys.argv
argc = len(param)

# columns 1 and 2: two nodes to form an event
# column 3: timestamp of the event

# SS = pd.read_csv('similarity_matrix_pop_to_pop_Np26.csv', header=0, index_col = 0)
# SS = pd.read_csv('similarity_matrix_pop_to_pop_Ng50_Np26.csv', header=0, index_col = 0)
SS = pd.read_csv('similarity_matrix_pop_to_pop_div_by_avg.csv', header=0, index_col = 0)
SS_values = SS.values
Ng = SS_values.shape[0] # number of genes
print(str(Ng) + " genes")

dg = pd.DataFrame(data=None, index=SS.index, columns=['x', 'y', 'd_from_ave']) # output, containing the result of MDS and the MU values
del SS
gc.collect()

print("The similarity matrix has been loaded.")

# print(SS.index)


Nf = SS_values.shape[1] # number of data columns for each gene
print(str(Nf) + " data columns") # There are 26 populations. So, there are 26*25/2 = 325 data columns

# dg = df.iloc[:, range(3,5)].copy() # deep copy
# df.drop(df.columns[range(7)], axis=1, inplace=True) # delete the first 7 columns and overwrite on df

if_use_Frobenius = 1

if if_use_Frobenius == 1:
    network_distance_matrix = np.zeros((Ng, Ng))
    for i in range(Ng):
    #    network_distance_matrix[i, i] = 0.0
        if i % 10 == 0:
            print(i)
        for j in range(i):
            network_distance_matrix[i, j] = Frobenius(SS_values[i,:], SS_values[j,:])
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
print("center = (", x_ave, y_ave, ")")
d_from_ave = np.sqrt(np.square(x - x_ave) + np.square(y - y_ave))
#df = df.drop(range(5, len(df.columns)), axis=1)
dg['x'] = x
dg['y'] = y
dg['d_from_ave'] = d_from_ave
# dg = dg.sort_values(by='d_from_ave', ascending=False)
# print(df)
dg.to_csv("MDS_result_alkan2.csv")

# pd.options.mode.chained_assignment = None  # default='warn'
