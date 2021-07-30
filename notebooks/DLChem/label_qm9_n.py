# -*- coding: utf-8 -*-

import fileinput
import os

import numpy as np

from schnetpack.datasets import QM9

from numpy import savetxt

class label_N():
    def __init__(self,qm9data,label_dir,number_data,element):
        self.label_dir = label_dir
        self.number_data = number_data
        self.qm9data = qm9data
        element = ' ' + element + '  '
        self.element = element

        label_n(self)
        
def label_n(self):

    element = self.element
    number_data = self.number_data
    label_dir = self.label_dir
    
    if not os.path.exists(label_dir):
        os.makedirs(label_dir)

    datan = np.zeros((1,2))    
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
                
                label = 0
                ##NOTE it is possible to have ammonium rarely TOTAL = 4
                if total == 3:
                    if countH == 3:
                        #ammonia NH3
                        label = 1
                    if countH == 2 and countC == 1:
                        #primary amine H2N-C
                        label = 2
                        
                        positions = connected_positions[j]
                        neighbor = min(positions)
                        index = positions.index(neighbor)
                        positions.pop(index)
                        
                        nn = neighboring_connections(name_mol,number_atoms,neighbor)
                        total2,countC2,countH2,countO2,countN2 = count_nn(nn,name_xyz)
                        

                        if total2 == 4: 
                            label = 3

                        
                        if total2 == 3:
                            label = 4
                            

                
                    if countH == 1 and countC == 2:
                        # secondary amine H2N-C2 
                        label = 5
                        
                        positions = connected_positions[j]
                        neighbor = min(positions)

                        index = positions.index(neighbor)
                        positions.pop(index)
                        
                        nn = neighboring_connections(name_mol,number_atoms,neighbor)
                        total2,countC2,countH2,countO2,countN2 = count_nn(nn,name_xyz)
                            
                        neighbor = min(positions)
                        
                        nn1 = neighboring_connections(name_mol,number_atoms,neighbor)
                        total3,countC3,countH3,countO3,countN3 = count_nn(nn1,name_xyz)

                        if total2 == 4 and total3 == 4:
                            label = 6
                        
                        if total2 == 3 or total3 == 3:
                            label = 7


                            
                    if countH == 1 and countN == 2:
                        label = 8
                    if countH ==1 and countN == 1 and countC ==1:
                        label = 9
                        ### check if neighbors neighbors is required by systematic graphing as well
                    if countC == 3:
                        label=10
                    if countC == 2 and countN == 1:
                        label = 11
                    if countC == 2 and countO == 1:
                        label = 12
                    if countC == 1 and countO == 1 and countN == 1:
                        label = 13
                    if countC == 1 and countN ==2:
                        label = 14
                if total == 2:
                    if countC == 2:
                        label = 15
                    if countC == 1 and countH == 1:
                        label = 16
                    if countC == 1 and countO == 1:
                        label = 17
                    if countC == 1 and countN == 1:
                        label= 18
                    if countO == 1 and countN == 1:
                        label = 19 
                    if countN == 2:
                        label = 20
                if total == 1:
                    if countC == 1:
                        label = 21
                if total == 4:
                    if countH == 3 and countC == 1:
                        label = 22
                    if countH == 2 and countC == 2:
                        label = 23

#                if label == 0:
#                    print('0!!!!!!')


                        
                row = np.array([[label,idx]])
                datan = np.vstack((datan,row))
                
    savetxt(label_dir+'label-n.csv',datan,delimiter=',')
    os.system('rm ' + label_dir + '*mol')
    os.system('rm ' + label_dir + '*xyz')
        
def write_coordinates_xyz_mol(self,idx):
    label_dir = self.label_dir
    qm9data = self.qm9data
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
                

