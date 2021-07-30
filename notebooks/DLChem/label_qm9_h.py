# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import fileinput
import os

import numpy as np

from schnetpack.datasets import QM9

from numpy import savetxt

import sys

import subprocess
import re

def get_platform():
    platforms = {
        'linux1' : 'Linux',
        'linux2' : 'Linux',
        'darwin' : 'OS X',
        'win32' : 'Windows'
    }
    if sys.platform not in platforms:
        return sys.platform
    
    return platforms[sys.platform]
    
class label_H():
    def __init__(self,qm9data,label_dir,number_data,element):
        self.label_dir = label_dir
        self.number_data = number_data
        self.qm9data = qm9data
        element = ' '+element + '  '
        self.element = element
        
        label_h(self)

        
def write_coordinates_xyz_mol(self,idx):
    qm9data=self.qm9data
    label_dir=self.label_dir

    name_xyz = label_dir + str(idx) + '.xyz'
    name_mol = label_dir + str(idx) + '.mol'
    
    at, props = qm9data.get_properties(idx)
    x = props['_positions'][:,0]
    y = props['_positions'][:,1]
    z = props['_positions'][:,2]
    x = x.numpy()
    y = y.numpy()
    z = z.numpy()

    file_xyz = open(name_xyz,mode='w',encoding='utf-8')

    number_atoms = len(z)
    number_atoms = int(number_atoms)
    file_xyz.write(str(len(z))+'\n')    
    file_xyz.write('title'+'\n')
    for i in range(number_atoms):
        if props['_atomic_numbers'][i] == 1:
            string = ('H' + ' ' +  str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
            file_xyz.write(string)
        if props['_atomic_numbers'][i] == 6:
            string = ('C' + ' ' +  str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
            file_xyz.write(string)
        if props['_atomic_numbers'][i] == 8:
            string = ('O' + ' ' +  str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
            file_xyz.write(string)
        if props['_atomic_numbers'][i] == 7:
            string = ('N' + ' ' +  str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
            file_xyz.write(string)
        if props['_atomic_numbers'][i] == 9:
            string = ('F' + ' ' +  str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
            file_xyz.write(string)
    file_xyz.close()
    os.system('obabel ' + name_xyz + ' -O ' + name_mol)
    return name_mol, props, number_atoms, name_xyz
    
def store_positions(self,idx,name_mol,props):
    element = self.element
    countmollines = 0
    array_of_element = []
    for line in fileinput.FileInput(name_mol,inplace=0):
        countmollines = countmollines + 1
        if element in line:
            array_of_element.append([])
            array_of_element[len(array_of_element)-1] = int(countmollines-4)
    return array_of_element


def connection_matrix(array_of_element,name_mol,number_atoms):
    count=0
    connected_positions = []
    array = []
    for j in range(len(array_of_element)):
        connected_positions.append([])
        array.append([])
        array[j] = str(array_of_element[j])
    for line in fileinput.FileInput(name_mol,inplace=0):
        count = count + 1
        if count >= number_atoms + 5:
            line = line.replace(line,line[:7])
            for j in range(len(array_of_element)):
                if ' ' + str(array_of_element[j]) + ' ' in line and 'RAD' not in line:
                    connected_positions[j].append([])
                    if str(line[1]+line[2]+line[3]) == ' '+array[j]+' ' or str(line[0]+line[1]+line[2]+line[3]) == ' '+array[j]+' ':
                        connected = str(line[4:6])
                        connected = int(connected)
                        connected_positions[j][len(connected_positions[j])-1]=connected
                    else:
                        connected = str(line[1]+line[2]+line[3])
                        connected = int(connected)
                        connected_positions[j][len(connected_positions[j])-1]=connected
#    print(connected_positions)
    return connected_positions


def neighboring_connections(name_mol,number_atoms,carbon_position):
    carbon_neighbor_connections = []
    countmollines = 0
    for line in fileinput.FileInput(name_mol,inplace=0):
        countmollines = countmollines + 1
        if countmollines >= number_atoms + 5:
            line = line.replace(line,line[:7])
            if ' ' + str(carbon_position) + ' ' in line and 'RAD' not in line:
                carbon_neighbor_connections.append([])
                if str(line[1]+line[2]+line[3]) == ' '+str(carbon_position)+' ' or str(line[0]+line[1]+line[2]+line[3]) == ' '+str(carbon_position)+' ':
                    connected = str(line[4:6])
                    connected = int(connected)
                    carbon_neighbor_connections[len(carbon_neighbor_connections)-1]=connected
                else:
                    connected = str(line[1]+line[2]+line[3])
                    connected = int(connected)
                    carbon_neighbor_connections[len(carbon_neighbor_connections)-1]=connected
    return carbon_neighbor_connections


def check_O(neighbor,name_xyz):
    countxyzlines = 0
    O_present=False

    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'O ' in line:
                O_present=True
    return O_present


def check_N(neighbor,name_xyz):
    countxyzlines = 0
    N_present=False

    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'N ' in line:
                N_present=True
    return N_present

def check_C(neighbor,name_xyz):
    countxyzlines = 0
    C_present=False

    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'C ' in line:
                C_present=True
    return C_present

def count_nn(nn,name_xyz):
    countC = 0
    countH = 0
    countO = 0
    countN = 0
    for j in range(len(nn)):
        target = nn[j]
        countxyzlines = 0
        for line2 in fileinput.FileInput(name_xyz,inplace=0):
            countxyzlines = countxyzlines+1
            if countxyzlines == nn[j]+2:
                if 'C ' in line2:
                    countC = countC+1
                if 'H ' in line2:
                    countH = countH+1
                if 'N ' in line2:
                    countN = countN+1
                if 'O ' in line2:
                    countO = countO+1
    total = countC + countH + countO + countN
    
    return total, countC, countH, countO, countN
                
        
def label_h(self):
    qm9data=self.qm9data
    label_dir=self.label_dir
    number_data=self.number_data
    element = self.element
    
    if not os.path.exists(label_dir):
        os.makedirs(label_dir)

    datah = np.zeros((1,2))
    for idx in range(number_data):
        load_bar = (idx/number_data)*(100)
        if load_bar%10 == 0:
            print(str(load_bar) +'% complete')
                
        # write xyz file and convert it to mol file
        name_mol, props, number_atoms, name_xyz = write_coordinates_xyz_mol(self,idx)
        
        if element == ' N  ':
            number = 7
        if element == ' O  ':
            number == 8
        if element == ' C  ':
            number = 6
        if element == ' H  ':
            number = 1
        if number in props['_atomic_numbers']:        
        #find and store all locations of nitrogens for this molecule
            positions = store_positions(self,idx,name_mol,props)
        
        
            # build bond matrix for each of these atoms
            connected_positions = connection_matrix(positions,name_mol,number_atoms)
            
            
            #BEGIN LABELLING 

            for j in range(len(connected_positions)):
                countC = 0
                countH = 0
                countO = 0
                countN = 0
                target = connected_positions[j]
                for i in range(len(target)):
                    countxyzlines=0
                    for line2 in fileinput.FileInput(name_xyz,inplace=0):
                        countxyzlines = countxyzlines+1
                        if countxyzlines == connected_positions[j][i]+2:
                            if 'C ' in line2:
                                countC = countC+1
                            if 'H ' in line2:
                                countH = countH+1
                            if 'N ' in line2:
                                countN = countN+1
                            if 'O ' in line2:
                                countO = countO+1
                total = countC + countH + countO + countN
                
                C_nbh = False
                N_nbh = False
                O_nbh = False
                
                if countC == 1:
                    C_nbh = True
                if countO == 1:
                    O_nbh = True
                if countN == 1:
                    N_nbh = True
                    
                positions = connected_positions[j]
                neighbor = min(positions)

                nn = neighboring_connections(name_mol,number_atoms,neighbor)
                total2,countC2,countH2,countO2,countN2 = count_nn(nn,name_xyz)
                
                label=0
                if C_nbh == True:
                    if total2 == 4:    
                        if countH2 == 4:
                            #this is methane, no cluster
                            label = 1
                        elif countH2 == 3 and countC2 == 1:
                            #these are all the methyl groups
                            label = 2
                        #methyl-O
                        elif countH2 == 3 and countO2 == 1:
                            label = 3
                        #methyl-N
                        elif countH2 == 3 and countN2 == 1:
                            label = 4
                        #this is methylene
                        elif countH2 == 2 and countC2 == 2:
                            label = 5
                        #methylene-N
                        elif countH2==2 and countC2==1 and countN2==1:
                            label = 6
                        elif countH2 == 2 and countC2 == 1 and countO2==1:
                            label = 7
                        elif countH2 == 2 and countO2 == 2:
                            label = 8
                        #methines
                        elif countH2==1 and countC2==3:
                            label=9
                        elif countH2==1 and countC2==2 and countN2==1:
                            label=10
                        elif countH2==1 and countC2==2 and countO2==1:
                            label=11
                        #outer acetal
                        elif countH2==1 and countC2==1 and countO2==2:
                            label=12 

                    elif total2 == 3:
                        #trigonal planar carbons
                        #arromatic
                        if countH2==1 and countC2==2:
                            label=13
                        elif countH2==1 and countC2==1 and countO2==1:
                            label=14
                        elif countH2==1 and countC2==1 and countN2==1:
                            label=15
                        elif countH2==1 and countO2==2:
                            label=16
                        elif countH2==1 and countN2==2:
                            label=17
                        elif countH2==1 and countN2==1 and countO2==1:
                            label=18
                        elif countH2==2 and countO2==1:
                            label=19
                            #formaldehyde
                    elif total2 == 2:
                        if countH2==1 and countC2==1:
                            label=20
                        elif countH2==1 and countN2==1:
                            label=21                     
                if N_nbh == True:
                    if total2 == 3:
                        #trigonal planar nitrogens
                        #first one is ammonia, no cluster
                        #ammonia            
                        if countH2==3:
                            label=22
                        #primary amine
                        elif countH2==2 and countC2==1:
                            
                            array_of_element = np.array([positions[0]])
                            connections = connection_matrix(array_of_element,name_mol,number_atoms)
                            positions2 = connections[0]
                            neighbor2 = min(positions2)

                            index = positions2.index(neighbor2)
                            positions2.pop(index)

                            nn2 = neighboring_connections(name_mol,number_atoms,neighbor2)
                            total3,countC3,countH3,countO3,countN3 = count_nn(nn2,name_xyz)

                            if total3 == 4:
                                label = 23

                            if total3 == 3:
                                label = 24

                        #secondary amine
                        elif countH2==1 and countC2==2:
                            label=25
                        elif countH2==1 and countN2==1 and countC2==1:
                            label=26
                        elif countH2==1 and countN2==2:
                            label=27
                    elif total2==2:
                        if countH2==1 and countC2==1:
                            label=28
                if O_nbh == True:
                    if total2==2:
                        #alcohols
                        if countH2==1 and countC2==1:

                            array_of_element = np.array([positions[0]])
                            connections = connection_matrix(array_of_element,name_mol,number_atoms)
                            positions2 = connections[0]
                            neighbor2 = min(positions2)

                            index = positions2.index(neighbor2)
                            positions2.pop(index)

                            nn2 = neighboring_connections(name_mol,number_atoms,neighbor2)
                            total3,countC3,countH3,countO3,countN3 = count_nn(nn2,name_xyz)

                            if total3 == 4:
                                label = 29

                            if total3 == 3:
                                label = 30

                        if countH2 == 2:
                            label=31
                        if countH2==1 and countN2==1:
                            label=32

                        
                row = np.array([[label,idx]])
                datah = np.vstack((datah,row))
                
    savetxt(label_dir+'label-h.csv',datah,delimiter=',')
    platform = get_platform()

    if platform == 'Linux':
        os.system('rm ' + label_dir + '*mol')
        os.system('rm ' + label_dir + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_dir = re.sub('../../../','/home/mobaxterm/Desktop/work/',label_dir)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*xyz"')

#qm9data = QM9('./qm9.db',download=True,remove_uncharacterized=True)
#label_dir = '/home/amerelsamman/Desktop/work/data/label/5000/'
#number_data = (5000) 
#element = ' H  '

#label_H(qm9data,label_dir,number_data,element)