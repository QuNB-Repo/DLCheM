# -*- coding: utf-8 -*-
"""
This subroutine studies changes in the representation of SchNet as it 
pertains to various molecular properties
"""

from schnetpack.datasets import QM9

import pandas as pd

import numpy as np
from numpy import genfromtxt

import matplotlib.pyplot as plt
import matplotlib

import math

def euclidean_diff(rep,number_data):
    diff = np.zeros((number_data-1,30))
    euc_dist = np.zeros((number_data-1))
    for idx in range(number_data-1):
    #random differences
        for i in range(30):
            diff[idx][i] = (rep[idx+1][i] - rep[idx][i])**2
        row = diff[idx]

        euc_dist[idx] = np.sum(row)
        euc_dist[idx] = math.sqrt(euc_dist[idx])
    return euc_dist

def plot2Dwlabel(x,y,labelx,labely,title,label,colors):
    fig1 = plt.figure(figsize=(5,4),dpi=100)
    ax1 = fig1.add_axes((0.1,0.1,0.9,0.9))
    ax1.set_xlabel(labelx)
    ax1.set_ylabel(labely)
    ax1.set_title(title)
    ax1.scatter(x,y,c=label['Target'],cmap=matplotlib.colors.ListedColormap(colors))

def plot2D(x,y,labelx,labely,title):
    fig1 = plt.figure(figsize=(5,4),dpi=100)
    ax1 = fig1.add_axes((0.1,0.1,0.9,0.9))
    ax1.set_xlabel(labelx)
    ax1.set_ylabel(labely)
    ax1.set_title(title)
    ax1.scatter(x,y)

def coefdif(rep,euc_const,coef,number_data):
    coefdif = np.zeros((number_data-1))
    
    for idx in range(number_data-1):
        if euc_const == False:
            coefdif[idx] = (rep[idx+1][coef] - rep[idx][coef])**2
#            coefdif[idx] = math.sqrt(coefdif[idx])
            
    return coefdif


rep = genfromtxt('../../data/diff/repfull_5000O.csv',delimiter=',')

number_data = 5000
euc_const = False             
euc_dist = euclidean_and_diff.euclidean_diff(rep,euc_const,number_data)

#the one for AE
coef = 30
AEdiff = euclidean_and_diff.coefdif(rep,euc_const,coef,number_data)
print(AEdiff)

#some other coefficient
coef = 29
coefdif = euclidean_and_diff.coefdif(rep,euc_const,coef,number_data)
print(coefdif)

labelx = 'Distance (30D)'
labely = '\u0394 AE (eV)'
title = 'Random Euclidean distances of schnet O representations in relation to AE change'

euclidean_and_diff.plot2D(euc_dist,AEdiff,labelx,labely,title)

coeff = 5
labelx = 'Coeff %s' %(coeff+1)
labely = 'AE (eV)' 
title = 'Coeff %s in relation to AE change' %(coeff+1)

label_file = '../../data/label_qm9/O/labelO5000.csv'
label = pd.read_csv(label_file,delimiter=',')
coef = rep[:,coeff]
AE = rep[:,30]
euclidean_and_diff.plot2Dwlabel(coef,AE,labelx,labely,title,label,colors)



#xlabel = '1st coeff'
#ylabel = 'AE (eV)'
#title = 'Coeff diff of schnet ether O representations in relation to AE diff'
#plot2D(coefdif,AEdiff,xlabel,ylabel,title)