# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 11:22:05 2021

@author: aelsamma
"""

import numpy as np
from numpy import genfromtxt, savetxt

import math 

global_map_data = '../../../data/global_map_data_network.csv'

colors = ['none','none','none','red','lightblue',
          'none','black','lightgreen','lime','none',
          'purple','coral','orange','none','none']

data = genfromtxt(global_map_data,delimiter=',')


# triangle matrix of distances between each pair of vectors
euc_dist_matrix = np.zeros((len(data),len(data)))
for i in range(len(data)):
    for j in range(len(data)):
        euc_dist_row = (data[i,:] - data[j,:])**2
        euc_dist_matrix[i][j] = math.sqrt(sum(euc_dist_row))


print(euc_dist_matrix)

euc_dist_matrix = euc_dist_matrix.flatten()

print(euc_dist_matrix)

count_first = 0
count_second = 0
network_data = np.zeros((len(data)**2,3))
for i in range(len(data)**2):
    network_data[i][0] = count_first
    network_data[i][1] = count_second
    network_data[i][2] = euc_dist_matrix[i]
    
    count_second = count_second+1
    if count_second == len(data):
        count_first = count_first+1
        count_second = 0
    
    
        
## prune the data, every len(data), delete 0 rows, then one, then 2, then 3,.... then 84 
count = 0
number_to_delete = 0 
done = False
number_deleted = 0 
for i in range(len(data)**2):
    if done == False: 
        for j in range(i,number_to_delete+i):
            network_data = np.delete(network_data,j-number_deleted-1,axis=0)
            number_deleted = number_deleted + 1
        done = True
    
    if count == len(data):
        count = 0
        number_to_delete = number_to_delete + 1
        done = False

    count = count + 1 
    
print(network_data)
savetxt('../../../data/network.csv',network_data,delimiter=',')
