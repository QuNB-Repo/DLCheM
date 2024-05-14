from input3 import*

import numpy as np
from numpy import genfromtxt, savetxt

#load dataset that will pick the transformation vectors
trans_data = genfromtxt(trans_filepath,delimiter=',',encoding='utf-8-sig')

#load init embedding vector dataset
init_data = genfromtxt(init_filepath,delimiter=',',encoding='utf-8-sig')


diff_vectors = np.zeros((len(trans_data),len(trans_data[0,:])))
for i in range(len(trans_data)):
    for j in range(len(trans_data[0,:])):
        diff_vectors[i][j] = trans_data[i][j] - init_data[i][j]

diff_vectors_mean = np.mean(diff_vectors,axis=0)

#add the vector mean to some target element init embedding what do you get in LDA?
target_data = genfromtxt(target_filepath, delimiter=',', encoding='utf-8-sig')

new_data = np.zeros((len(target_data),len(target_data[0,:])))
for i in range(len(target_data)):
    new_data[i,:] = target_data[i,:] + diff_vectors_mean

print(new_data)

savetxt('frominit.csv',new_data,delimiter=',')
