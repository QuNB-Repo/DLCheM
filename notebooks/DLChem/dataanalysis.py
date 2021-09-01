# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 18:02:03 2021

@author: aelsamma
"""

import numpy as np
from numpy import genfromtxt, savetxt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import  PCA
from sklearn.manifold import TSNE


import scipy.linalg as la

import statistics
import math

def pca(data_file_path,n_components,scale_data):


    data = genfromtxt(data_file_path +'.csv',delimiter=',')
    
    if scale_data == True:
        scaler = StandardScaler()
        scaler.fit(data)
        data = scaler.transform(data)
    #perform PCA decomposition of the data
        savetxt(data_file_path + 'scaled.csv', data, delimiter=',')
        
    
    pca = PCA(n_components)
    pca.fit(data)
    x_pca = pca.transform(data)
    
    cov = pca.get_covariance()
    eig, ev = la.eig(cov)
    total = sum(eig.real)
    evt=np.transpose(ev)
    unit=np.matmul(ev,evt)
    #plot the PCA if dimension 2 or 3, if element H or O for colors/legend
    pca_file = data_file_path +'pca.csv'
    savetxt(pca_file, x_pca, delimiter=',')

    
    return x_pca, eig, ev, cov

    
def tsne(data_file_path,dimension,perp):
    
    data = genfromtxt(data_file_path+'.csv',delimiter=',')
    
    X = TSNE(n_components=dimension,perplexity=perp).fit_transform(data)

    tsne_file = data_file_path + 'tsne.csv'
    savetxt(tsne_file, X, delimiter=',')
    
    return X          
    

def pearson(data_file,ae_data_file):
    ae_data = genfromtxt(ae_data_file, delimiter=',')
    data = genfromtxt(data_file,delimiter=',')
   
    mean_data = np.zeros((30))
    stdev_data = np.zeros((30))
    for i in range((30)):
        mean_data[i] = statistics.mean(data[:,i])
        stdev_data[i] = statistics.stdev(data[:,i])
        
    
    
    mean_ae = statistics.mean(ae_data[:,2])
    stdev_ae = statistics.stdev(ae_data[:,2])
    
    pear = np.zeros((30))
    covar = np.zeros((30))
    for i in range((30)):
        for j in range(len(data)):
                covar[i] = covar[i] + (data[j][i]-mean_data[i])*(ae_data[j][2]-mean_ae)
    
    covar = covar/len(data)
    
    for i in range((30)):
        pear[i] = covar[i]/(stdev_data[i]*stdev_ae)
    
    return pear
    
