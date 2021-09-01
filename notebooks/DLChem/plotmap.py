# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 09:53:05 2021

@author: aelsamma
"""

import pandas as pd
from schnetpack.datasets import QM9

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from numpy import genfromtxt

from mpl_toolkits import mplot3d


qm9data = QM9('../../../data/QM9/qm9.db', download=False,remove_uncharacterized=True)
global_map_data = '../../../data/global_map_data.csv'

colors = ['none','none','none','red','lightblue',
          'none','black','lightgreen','lime','none',
          'purple','coral','orange','none','none']

data = genfromtxt(global_map_data,delimiter=',')

fig = plt.figure(figsize=(5,4),dpi=100)
ax = plt.axes(projection='3d')
ax.scatter(data[:,0],data[:,1],c=data[:,30],cmap=matplotlib.colors.ListedColormap(colors))

ax.view_init(30,290)

calc6561 = np.zeros((3))
for i in range(3):
    calc6561[i] = data[60][i] - data[64][i]
    
print(calc6561)

add = np.zeros((85,1))
for j in range(85):
    for i in range(4,5):
        add[j] = add[j] + data[j][i]
    

print(add)
print(len(add))
print(len(data[:,31]))

fig2 = plt.figure(figsize=(5,4),dpi=100)
ax2 = plt.axes()
ax2.scatter(add,data[:,31],c=data[:,30],cmap=matplotlib.colors.ListedColormap(colors))


column = data[:,4]
column = np.sort(column, axis=-1)
print(column)


