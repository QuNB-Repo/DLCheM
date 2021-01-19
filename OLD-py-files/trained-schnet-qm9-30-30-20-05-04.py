#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loads the best energy model after schnet training.
user evaluates the model on the QM9 database



@author: amerelsamman
"""


import schnetpack as spk
import torch
import schnetpack.nn 
import schnetpack.data
from schnetpack.datasets import QM9

# Load best model where it was stored during the training

model = torch.load("./qm9tut/best_model", map_location=torch.device('cpu'))

qm9data = QM9('./qm9.db', download=True, remove_uncharacterized=True)

train, val, test = spk.data.train_test_split(qm9data, split_file='./qm9tut/split.npz')

device = 'cpu'

converter = spk.data.AtomsConverter(device=device)

from xlwt import Workbook


#Control excel sheets to output data

wb1 = Workbook()
sheet1 = wb1.add_sheet('Sheet1')

#Write the number of molecules you will evaluate
sheet1.write(0,0, 100)
sheet1.write(0,1, 40)

#Some important counters for coding output into excel sheets 
index=0


#BEGIN AUTOMATED CALCULATIONS, choose the range of QM9 to perform the calculations in
for idx in range(1334,1335):

    test_loader = spk.AtomsLoader(test, batch_size=100)
    converter = spk.data.AtomsConverter(device=device)
    at, props = qm9data.get_properties(idx)

    x = props['_positions'][ :,0]
    y = props['_positions'][ :,1]
    z = props['_positions'][ :,2]


    inputs = converter(at)
    print('Keys:', list(inputs.keys()))
    print('Truth:', props[QM9.U0].cpu().numpy()[0])

    calculator = spk.interfaces.SpkCalculator(model=model, device=device, energy=QM9.U0)
    at.set_calculator(calculator)


    print('Prediction:', at.get_total_energy())
    print('idx',idx)


    import numpy as np
 

#Output atomization energies from schnetpack (you have to first set the variable in schnetpack as global)
    from schnetpack.atomistic import output_modules

    atomenergies = (output_modules.yi)

    atomenergies = atomenergies.detach().numpy()
    atomenergies = atomenergies.astype(float)
    print(atomenergies[0][1])

    import math 

    CMat = np.zeros((len(z),len(z)))

    Z = props['_atomic_numbers']

    for i in range(len(z)):
        for j in range(len(z)):
            CMat[i,j] = CMat[i,j] + math.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)
            if i == j:
                CMat[i,j]=0.5*Z[i]**(2.4)
            if i != j:
                CMat[i,j]=Z[i]*Z[j]/CMat[i,j]
            
    Z=Z.cpu().detach().numpy() 
    import scipy.linalg as la
 
    eig, ev = la.eig(CMat)
    eigen = eig.astype(float)
    
    
    #extract atomenergies into a 1D array
    AE = np.zeros(len(z))
    for i in range(len(z)):
        AE[i] = atomenergies[0][i]
    totz = 0
    for i in range(len(z)):
        totz = totz + Z[i]
        totz = float(totz)
    print('totz',totz)
    print('atomic numbers', Z)
    #write atomization energies of the first 20 atoms, if there is no atom, just write 0
 #   if len(z) == 16:
 #       if totz == 60
    if len(z) == 9:
        if totz == 50:
            sheet1.write(index+1,0,idx)
            count = 0 
            for i in range(len(z)):
                if Z[i] == 6:
                    sheet1.write(index+1,count+1, AE[i])
                    count=count+1
            for i in range(len(z)):
                if Z[i] == 1:
                    sheet1.write(index+1,count+1,AE[i])
                    count=count+1
            for i in range(len(z)):
                if Z[i] == 7:
                    sheet1.write(index+1,count+1, AE[i])
                    count=count+1
            for i in range(len(z)):
                if Z[i] == 8:
                    sheet1.write(index+1,count+1, AE[i])
                    count=count+1
            for i in range(len(z)):
                if Z[i] == 9:
                    sheet1.write(index+1,count+1, AE[i])
                    count=count+1
            for i in range(len(z)):
                sheet1.write(index+1,len(z)+i+1,eigen[i])
            index=index+1
    wb1.save('./pca/ALL-5000.xls') 
    x = props['_positions'][ :,0]
    y = props['_positions'][ :,1]
    z = props['_positions'][ :,2]
    x = x.numpy()
    y = y.numpy()
    z = z.numpy()

    j = len(props['_positions'])
    print('number of atoms', j)
    for i in range(j):
        if props['_atomic_numbers'][i] == 1:
            print('H',x[i],y[i],z[i])
        if props['_atomic_numbers'][i] == 6:
            print('C',x[i],y[i],z[i])
        if props['_atomic_numbers'][i] == 7:
            print('N',x[i],y[i],z[i])     
        if props['_atomic_numbers'][i] == 8:
            print('O',x[i],y[i],z[i])
        if props['_atomic_numbers'][i] == 9:
            print('F',x[i],y[i],z[i]) 









#for computing the average atomization of each element in a molecule
#    numH = 0
#    numC = 0
#    numN = 0
#    numO = 0
#    numF=0
#    h=0
#    c=0
#    n=0
#    o=0
#    f=0
#    for i in range(len(z)):
#        atomicnumbers = props['_atomic_numbers'][i]
#        if atomicnumbers == 1:
#            numH = numH + atomenergies[0,i]
#            h = h +1
#        if atomicnumbers == 6:
#            numC = numC + atomenergies[0,i]
#            c = c +1
#        if atomicnumbers == 7:
#            numN = numN + atomenergies[0,i]
#            n = n+1
#        if atomicnumbers == 8:
#            numO = numO + atomenergies[0,i]
#            o = o + 1
#        if atomicnumbers == 9:
#            numF = numF + atomenergies[0,i]
#            f = f + 1
#    if numH != 0:
#        aveH = numH/h
#    else:
#        aveH = -13.61312
#    if numC != 0:
#        aveC = numC/c
#    else:
#        aveC = -1029.86302
#    if numN != 0:
#        aveN = numN/n
#    else:
#        aveN = -1485.302400 
#    if numO != 0:
#        aveO = numO/o
#    else:
#        aveO = -2042.6110546
#    if numF !=0:
#        aveF = numF/f
#    else:
#        aveF = -2713.48462395
#    print(aveC)
#    print(aveH)
#    print(aveN)
#    print(aveO)
#    print(aveF)
#    aveH = float(aveH)
#    aveC = float(aveC)
#    aveN = float(aveN)
#    aveO = float(aveO)
 #   aveF = float(aveF)
#    if len(z) == 14:
#        sheet1.write(index+1,0, aveH)
#        sheet1.write(index+1,1, aveC)
#        sheet1.write(index+1,2, aveN)
#        sheet1.write(index+1,3, aveO)
#        sheet1.write(index+1,4, aveF)
#        for i in range(len(z)):
#            sheet1.write(index+1,5+i, eigen[i])
#        index=index+1
    
#    import matplotlib.pyplot as plt
#    import matplotlib.


# Data for three-dimensional scattered points
#    for i in range(len(z)):
#        zdata = z[i]
#        ydata = y[i]
#        xdata = x[i]
#        atomicnumbers = props['_atomic_numbers'][i]
#        if atomicnumbers == 6:
#            c = 'g'
#        if atomicnumbers == 1:
#            c = 'b'
#        if atomicnumbers == 7:
#            c = 'y'
#        if atomicnumbers == 8:
#            c = 'r'
#        if atomicnumbers == 9:
#            c ='m'
#        ax.scatter3D(xdata, ydata, zdata, c=c)
#    ax.set_xlabel('X Label')
#    ax.set_ylabel('Y Label')
#    ax.set_zlabel('Z Label')
#    plt.show()