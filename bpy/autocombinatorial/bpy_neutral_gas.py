#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 05:41:00 2021

@author: amerelsamman
"""


import numpy as np
import fileinput
import random

num = []
for i in range(20):
    num.append(i)

print(num)
for i in range(20):
    bpy = 'bpy_neutral_gas.inp'
    string1 = ('H     1.0     3.35336     3.40622     0.60958',
               'Cl    17.0     3.90177     3.09707     0.62015',
               'C     6.0     3.84735     3.30245     0.27476' + '\n' 'H     1.0     3.68110     2.22497     0.16954' + '\n' 'H     1.0     4.32605     3.47028     1.24547' + '\n' 'H     1.0     4.54404     3.60780    -0.51341',
               'N     7.0     3.86509     3.03339     0.05557' + '\n' 'C     6.0     3.76803     1.57378    -0.08430' + '\n' 'H     1.0     2.95536     1.16318     0.52427' + '\n' 'H     1.0     4.68489     1.08022     0.25695' + '\n' 'H     1.0     3.60577     1.30298    -1.13331' + '\n' 'C     6.0     5.22551     3.52073    -0.19046' + '\n' 'H     1.0     5.37533     4.57331     0.05407' + '\n' 'H     1.0     5.49047     3.37600    -1.24326' + '\n' 'H     1.0     5.94155     2.96168     0.42322',
               'C     6.0     3.91501     3.07784    -0.12209' + '\n' 'N     7.0     4.86353     2.41388    -0.19908',
               'N     7.0     3.84558     3.23902     0.23632' + '\n' 'O     8.0     3.76638     1.95281    -0.28324' + '\n' 'O     8.0     4.98996     3.94911    -0.10550',
               'C     6.0     4.19686     3.08319    -0.21833' + '\n' 'C     6.0     4.19701     1.60092    -0.48964' + '\n' 'H     1.0     3.96570     1.05940     0.43051' + '\n' 'H     1.0     5.19438     1.30539    -0.83275' + '\n' 'H     1.0     3.49652     1.33893    -1.28492' + '\n' 'O     8.0     5.26759     3.67531    -0.05997',
               'C     6.0     4.24118     2.97266    -0.35506' + '\n' 'O     8.0     4.16501     1.66896    -0.01123' + '\n' 'H     1.0     3.37076     1.43930     0.50280' + '\n' 'O     8.0     5.27475     3.48653    -0.75289',
               'C     6.0     4.36368     2.73940    -0.40260' + '\n' 'O     8.0     4.13541     1.41840    -0.62602' + '\n' 'O     8.0     5.47311     3.25417    -0.33560' + '\n' 'C     6.0     5.33123     0.65471    -0.79673' + '\n' 'H     1.0     5.05170    -0.39492    -0.92086' + '\n' 'H     1.0     5.97427     0.73819     0.08563' + '\n' 'H     1.0     5.86531     0.97967    -1.69562')
    string2 = ('H     1.0    -2.69741     0.60890    -0.60830',
               'Cl    17.0    -2.73807    -0.01290    -0.74258',
               'C     6.0    -2.64738     0.15360    -0.64525' + '\n' 'H     1.0    -3.59123    -0.07620    -1.15121' + '\n' 'H     1.0    -2.62401    -0.40152     0.29838' + '\n' 'H     1.0    -1.84088    -0.21877    -1.28638',
               'N     7.0    -2.60687     0.24571    -0.54764' +'\n' 'C     6.0    -3.90476    -0.42097    -0.45612' +'\n' 'H     1.0    -4.60495     0.11655     0.19095' + '\n' 'H     1.0    -3.79455    -1.42301    -0.02515' + '\n' 'H     1.0    -4.34446    -0.52337    -1.45339' + '\n' 'C     6.0    -1.49993    -0.66783    -0.84020' + '\n' 'H     1.0    -0.67147    -0.16675    -1.35215' + '\n' 'H     1.0    -1.82833    -1.47519    -1.50468' + '\n' 'H     1.0    -1.12562    -1.11350     0.08701',
               'C     6.0    -2.51276     0.17881    -0.58088' + '\n' 'N     7.0    -2.62868    -0.96833    -0.71277',
               'N     7.0    -2.59940     0.15819    -0.79594' + '\n' 'O     8.0    -3.83162    -0.37090    -0.43197' + '\n' 'O     8.0    -1.49767    -0.62516    -0.47466',
               'C     6.0    -2.45698    -0.07601    -0.63992' + '\n' 'C     6.0    -1.27160    -0.98908    -0.82119' + '\n' 'H     1.0    -0.77520    -0.76964    -1.76945' + '\n' 'H     1.0    -1.62244    -2.02618    -0.84647' + '\n' 'H     1.0    -0.57801    -0.89716     0.01742' + '\n' 'O     8.0    -3.59875    -0.53775    -0.65855',
               'C     6.0    -2.23889    -0.20194    -0.64483' + '\n' 'O     8.0    -1.22879    -0.79910    -1.31407' + '\n' 'H     1.0    -0.65056    -0.17731    -1.79141' + '\n' 'O     8.0    -3.20311    -0.81784    -0.22029',
               'C     6.0    -2.21214    -0.35444    -0.53693' + '\n' 'O     8.0    -1.06365    -1.03375    -0.27785' + '\n' 'O     8.0    -3.28701    -0.87250    -0.80849' + '\n' 'C     6.0    -1.20323    -2.45427    -0.32786' + '\n' 'H     1.0    -0.28380    -2.89824     0.06363' + '\n' 'H     1.0    -1.33554    -2.78111    -1.36358' + '\n' 'H     1.0    -2.03766    -2.79272     0.29521')

    new_file=open(num[i],mode="w",encoding="utf-8")
    for line in fileinput.FileInput(bpy,inplace=0):
        if 'H     1.0     3.35336     3.40622     0.60958' in line:
            choice = random.choice(string1)
            line = line.rstrip()
            line = line.replace(line, choice +'\n')
        if 'H     1.0    -2.69741     0.60890    -0.60830' in line:
            choice = random.choice(string2)
            line = line.rstrip()
            line = line.replace(line, choice + '\n')
        print(line, end='')
        new_file.write(line)
    
new_file.close()