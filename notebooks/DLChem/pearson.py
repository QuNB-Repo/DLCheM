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
    
    print(data)
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
    
    stdev = statistics.stdev(data[1])

    #calculate pearson coefficients
    pearson = np.zeros((len(cov),len(cov)))
    for i in range(len(cov)):
        for j in range(len(cov)):
            stdevx = statistics.stdev(data[i])
            stdevy = statistics.stdev(data[j])
            pearson[i][j]=cov[i][j]/(stdevx*stdevy)
    
    print(pearson)    
    savetxt('../../../data/pca/'+name_data+'pca.csv',x_pca,delimiter=',')
#    pearson_sums = np.zeros((30))
#    for i in range(30):
#        squ = pearson[j][0]**2
#        for j in range(1,30):
#            pearson_sums[i] = squ + pearson[i][j]**2
    
#    print(pearson_sums)
    
element = 'O'
name_data='rep-5000-model1H-quint'
number_data=5000
file_path = '../../../data/schnet/' + name_data + '.csv' 
label_path = '../../../data/label_qm9/%s/label%s5000.csv' %(element,element) 
scale_data = True
n_components=30
dimension=2
pearson(element,file_path,label_path,scale_data,n_components,dimension)    


