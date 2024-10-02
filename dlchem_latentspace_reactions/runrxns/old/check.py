from checkinput import *
import numpy as np
from numpy import genfromtxt, savetxt

# concatenate the embedding vectors of fgtransform file with the embedding vectors of extractembeddings file for an overall comparison
true_vecs = genfromtxt(truetransformed_filepath,delimiter=',',encoding='utf-8-sig')
extract_vecs = genfromtxt(extractembs_filepath,delimiter=',',encoding='utf-8-sig')
pred_vecs = genfromtxt(predtransformed_filepath,delimiter=',',encoding='utf-8-sig')


full_vecs = np.vstack((true_vecs,extract_vecs))


print(len(true_vecs))
# take norm difference between first embedding and all the rest and find out which is closest 

lowest_distance = 100
index = np.zeros((len(pred_vecs)))
for j in range(len(pred_vecs)):
    for i in range(len(true_vecs)):
        distance = np.linalg.norm(pred_vecs[j,:]-true_vecs[i,:])
#    print(distance)
        if distance < lowest_distance:
            lowest_distance = distance 
            ind = i
            save_vec = true_vecs[i,:]

    index[j] = ind

print(index)