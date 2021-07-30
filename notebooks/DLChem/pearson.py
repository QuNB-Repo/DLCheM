# -*- coding: utf-8 -*-

import numpy as np
from numpy import genfromtxt, savetxt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import  PCA


import matplotlib.pyplot as plt
import matplotlib


import scipy.linalg as la

import statistics


def pearson(element,file_path,label_path,scale_data,n_components,dimension):
    if element == 'H':
         colors=['purple','none','palegreen','pink','lightblue',
             'mediumseagreen','lightskyblue','lightcoral','none','darkseagreen',
             'cornflowerblue','indianred','none','yellowgreen','palevioletred',
             'aquamarine','none','steelblue','mediumvioletred','none',
             'darkolivegreen','none','none','khaki','hotpink',
             'navajowhite','none','none','thistle','darksalmon',
             'tomato','none','burlywood']
    if element == 'C':
        colors=['black','none','#ccffcc','#ccffff','#ffcccc',
            '#99ff66','none','none','#66ffff','#ff9999',
            '#77b300','#33cccc','#ff6666','none','#666633',
            '#cc0000','none','none','none','none',
            'none','#ffffcc','#ffcc99','#66ccff','#ccccff',
            'none','#00cc66','#009999','#ff9933','#0099ff',
            '#9999ff','black','#00ff99','#006666','black',
            'none','#99ccff', '#3333cc', '#ff6600','none',
            'blue','none','none','black','black',
            'black','none']
    if element == 'N':
        colors=['none','none','none','palegreen','khaki',
                'none','mediumaquamarine','mediumseagreen','none','orange',
                'darkseagreen','none','none','none','none',
                'lightskyblue', 'steelblue','lightcoral','none','none',
                'tomato','plum','purple','hotpink']
    
    if element == 'O':
        colors=['none','none','tomato','none','mediumaquamarine',
               'mediumseagreen','palegreen','none','lightblue','none',
                'lightskyblue','none','cornflowerblue','khaki','NONE',
                'none','lightpink','mistyrose','deeppink','hotpink',
                'mediumvioletred','plum','purple','NONE','mediumorchid']

    data = genfromtxt(file_path,delimiter=',')
    
    if scale_data == True:
        scaler = StandardScaler()
        scaler.fit(data)
        data = scaler.transform(data)
    #perform PCA decomposition of the data

    
    pca = PCA(n_components)
    pca.fit(data)
    x_pca = pca.transform(data)
    
    cov = pca.get_covariance()


    eig, ev = la.eig(cov)
    eig = eig.real 
    print(eig)
    
    stdev = statistics.stdev(data[1])

    #calculate pearson coefficients
    pearson = np.zeros((len(cov),len(cov)))
    for i in range(len(cov)):
        for j in range(len(cov)):
            stdevx = statistics.stdev(data[i])
            stdevy = statistics.stdev(data[j])
            pearson[i][j]=cov[i][j]/(stdevx*stdevy)
        

    print(pearson[0])
    
    pca_file = '../../../data/pca/%s.csv' %(str(number_data)+name_data + element)
    savetxt(pca_file, x_pca, delimiter=',')
    
    #-1 from pca file
    #propanol-ethanol
    diff1 = x_pca[173]-x_pca[51]
    
    #ethanol - methanol
    diff2 = x_pca[51]-x_pca[23]
    
    index1 = np.argmax(diff1)
    index2 = np.argmax(diff2)
    print(index1)
    print(index2)

    
    
element = 'H'
name_data='int0'
number_data=5000
file_path = '../../../data/schnet/5000/' + name_data + '%s.csv' %(element)
label_path = '../../../data/label_qm9/%s/5000/label%s.csv' %(element,element) 
scale_data = True
n_components=30
dimension=2
pearson(element,file_path,label_path,scale_data,n_components,dimension)    


