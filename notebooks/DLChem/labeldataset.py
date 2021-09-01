# -*- coding: utf-8 -*-


import fileinput
import os

import numpy as np

from schnetpack.datasets import QM9

from numpy import savetxt

import sys
import subprocess
import re

from DLChem.utils import utils


class label():
    def __init__(self,qm9data,label_dir,label_name,number_data,element):
        self.label_dir = label_dir
        self.number_data = number_data
        self.qm9data = qm9data
        element = ' ' + element + '  '
        self.element = element
        self.label_name = label_name
        
        if element == ' O  ':
            label_o(self)
        if element == ' N  ':
            label_n(self)
        if element == ' H  ':
            label_h(self)
        
        
def label_o(self):
    number_data = self.number_data
    element = self.element
    label_dir = self.label_dir
    label_name = self.label_name
    qm9data = self.qm9data
    
    if not os.path.exists(label_dir):
        os.makedirs(label_dir)
    
    datao = np.zeros((1,2))
    for idx in range(number_data):
        load_bar = (idx/number_data)*(100)
        if load_bar%10 == 0:
            print(str(load_bar) +'% complete')
            
        # write xyz file and convert it to mol file
        at, props = qm9data.get_properties(idx)
        number_atoms = len(props['_atomic_numbers'])
        
        name_xyz = utils.write_xyz_from_db(props,idx,label_dir)
        name_mol = utils.xyz_to_mol(idx,label_dir)
        
        
        if 8 in props['_atomic_numbers']:
            
        #find and store all locations of nitrogens for this molecule
            positions = utils.store_positions(idx,name_mol,props,element)
        
        
            # build bond matrix for each of these atoms
            connected_positions = utils.connection_matrix(positions,name_mol,number_atoms)
            
            
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
                

                
                if total == 2:
                    if countH == 2:
                        label = 1
                    
                    if countC == 2:
                        label = 2
                    
                    if countC == 1 and countN == 1:
                        label = 3

                    if countH == 1 and countC == 1:
                        
                        positions = connected_positions[j]
                        neighbor = min(positions)
                        index = positions.index(neighbor)
                        positions.pop(index)
                        
                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
                        total2,countC2,countH2,countO2,countN2 = utils.count_nn(nn,name_xyz)
                        

                        if total2 == 4: 
                            if countC2 == 3:
                                label = 4
                            if countC2 == 2 and countH2 == 1:
                                label = 5
                            if countC2 == 1 and countH2 == 2:
                                label = 6
                            if countH2 == 3:
                                label = 7
                      
                        if total2 == 3:

                            if countC2 == 2:
                                label = 8
                            if countC2 == 1 and countO2 == 2:
                                label = 9
                            if countC2 == 1 and countN2 == 1:
                                label = 10
                            if countO2 == 2 and countN2 == 1:
                                label = 11
                            if countN2 == 2:
                                label = 12

                         
                    if countH == 1 and countN == 1:
                        label = 13
                    if countN == 2:
                        label = 14
                if total == 1:
                    if countC == 1:
                        positions = connected_positions[j]
                        neighbor = min(positions)
                        index = positions.index(neighbor)
                        positions.pop(index)
                        
                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
                        total2,countC2,countH2,countO2,countN2 = utils.count_nn(nn,name_xyz)
                        
                        if total2 == 3: 
                            if countH2 == 2:
                                label = 15
                            if countC2 == 2:
                                label = 16
                            if countC2 == 1 and countH2 == 1:
                                label = 17
                            if countC2 == 1 and countN2 == 1:
                                label = 18
                            if countC2 == 1 and countO2 == 2:
                                label = 19
                            if countN2 == 1 and countH2 == 1:
                                label = 20
                            if countN2 == 1 and countO2 == 2:
                                label = 21
                            if countN2 == 2:
                                label = 22
                            if countO2 == 3:
                                label = 23
                            if countO2 == 2 and countH2 == 1:
                                label = 24
                            
                if label == 0:
                    print('hereZERO!!!!')

                        
                row = np.array([[label,idx]])
                datao = np.vstack((datao,row))
                
    savetxt(label_dir+label_name,datao,delimiter=',')


    platform = utils.get_platform()
    
    #rtied a million times not to use MOBA... there is no combatibility. MOBA must be used at the end because it needs to be closed
    #has to do with the backslash in windows vs how it is used here
    if platform == 'Linux':
        os.system('rm ' + label_dir + '*mol')
        os.system('rm ' + label_dir + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_dir = re.sub('../../','/home/mobaxterm/Desktop/DLChem/',label_dir)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*xyz"')

def label_n(self):

    element = self.element
    number_data = self.number_data
    label_dir = self.label_dir
    label_name = self.label_name
    qm9data = self.qm9data
    
    if not os.path.exists(label_dir):
        os.makedirs(label_dir)

    datan = np.zeros((1,2))    
    for idx in range(number_data):
        load_bar = (idx/number_data)*(100)
        if load_bar%10 == 0:
            print(str(load_bar) +'% complete')
        
        # write xyz file and convert it to mol file
        at, props = qm9data.get_properties(idx)
        number_atoms = len(props['_atomic_numbers'])
        
        name_xyz = utils.write_xyz_from_db(props,idx,label_dir)
        name_mol = utils.xyz_to_mol(idx,label_dir)
        
        
        if 7 in props['_atomic_numbers']:        
        #find and store all locations of nitrogens for this molecule
            positions = utils.store_positions(idx,name_mol,props,element)
        
            # build bond matrix for each of these atoms
            connected_positions = utils.connection_matrix(positions,name_mol,number_atoms)
            
            
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
                        
                        positions = connected_positions[j]
                        neighbor = min(positions)
                        index = positions.index(neighbor)
                        positions.pop(index)
                        
                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
                        total2,countC2,countH2,countO2,countN2 = utils.count_nn(nn,name_xyz)
                        

                        if total2 == 4: 
                            if countC2 == 3:
                                label = 2
                            if countC2 == 2 and countO2 == 1:
                                label = 3
                            if countC2 == 2 and countN2 == 2:
                                label = 4
                            if countC2 == 1 and countO2 == 2:
                                label = 5
                            if countC2 == 1 and countN2 == 3:
                                label = 6
                            if countC2 == 1 and countO2 == 1 and countN2 == 2:
                                label = 7
                            if countO2 == 3:
                                label = 8
                            if countN2 == 4:
                                label = 9
                            if countN2 == 3 and countO2 == 1:
                                label = 10
                            if countN2 == 2 and countO2 == 2:
                                label = 11
                            if countH2 == 1 and countC2 == 2:
                                label = 12
                            if countH2 == 1 and countC2 == 1 and countN2 == 1:
                                label = 13
                            if countH2 == 1 and countC2 == 1 and countO2 == 1:
                                label = 14
                            if countH2 == 1 and countN2 == 2:
                                label = 15
                            if countH2 == 1 and countO2 == 2:
                                label = 16
                            if countH2 == 1 and countO2 == 1 and countN2 == 1:
                                label = 17
                            if countH2 == 2 and countC2 == 1:
                                label = 18
                                print('18!!!!!!!!!!!')
                                print('youforgotNneeddifftreatment,mustaddtototal2')
                                print(idx)
                            if countH2 == 2 and countO2 == 1:
                                label = 19
                            if countH2 == 2 and countN2 == 1:
                                label = 20
                            if countH2 == 3:
                                label = 21
                                                            
                        if total2 == 3:
                            if countC2 == 2:
                                label = 22
                            if countC2 == 1 and countO2 == 1:
                                label = 23
                            if countO2 == 2:
                                label = 24
                            if countO2 == 1 and countN2 == 2:
                                label = 25
                            if countN2 == 3: 
                                label = 26
                            if countN2 == 2 and countC2 == 1:
                                label = 27
                            if countC2 == 1 and countH2 == 1:
                                label = 28
                            if countN2 == 2 and countH2 == 1:
                                label = 29
                            if countO2 == 1 and countH2 ==1:
                                label = 30
                                
                                

                
                    if countH == 1 and countC == 2:
                        # secondary amine H2N-C2 
                        label = 31
                        
#                        positions = connected_positions[j]
#                        neighbor = min(positions)

#                        index = positions.index(neighbor)
#                        positions.pop(index)
                        
#                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
#                        total2,countC2,countH2,countO2,countN2 = utils.count_nn(nn,name_xyz)
                            
#                        neighbor = min(positions)
                        
#                        nn1 = utils.neighboring_connections(name_mol,number_atoms,neighbor)
#                        total3,countC3,countH3,countO3,countN3 = utils.count_nn(nn1,name_xyz)

#                        if total2 == 4 and total3 == 4:
#                            label = 24
                        
#                        if total2 == 3 or total3 == 3:
#                            label = 25


                            
                    if countH == 1 and countN == 2:
                        label = 32
                    if countH ==1 and countN == 1 and countC ==1:
                        label = 33
                        ### check if neighbors neighbors is required by systematic graphing as well
                    if countC == 3:
                        label=34
                    if countC == 2 and countN == 1:
                        label = 35
                    if countC == 2 and countO == 1:
                        label = 36
                    if countC == 1 and countO == 1 and countN == 1:
                        label = 37
                    if countC == 1 and countN ==2:
                        label = 38
                if total == 2:
                    if countC == 2:
                        label = 39
                    if countC == 1 and countH == 1:
                        label = 40
                    if countC == 1 and countO == 1:
                        label = 41
                    if countC == 1 and countN == 1:
                        label= 42
                    if countO == 1 and countN == 1:
                        label = 43 
                    if countN == 2:
                        label = 44
                if total == 1:
                    if countC == 1:
                        label = 45
                if total == 4:
                    if countH == 3 and countC == 1:
                        label = 46
                    if countH == 2 and countC == 2:
                        label = 47

                if label == 0:
                    print('ZERO!!!!!!')
                    print(idx)


                        
                row = np.array([[label,idx]])
                datan = np.vstack((datan,row))
                
    savetxt(label_dir+label_name,datan,delimiter=',')

    platform = utils.get_platform()
    
    if platform == 'Linux':
        os.system('rm ' + label_dir + '*mol')
        os.system('rm ' + label_dir + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_dir = re.sub('../../','/home/mobaxterm/Desktop/DLChem/',label_dir)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*xyz"')
        
        
def label_h(self):
    qm9data=self.qm9data
    label_dir=self.label_dir
    number_data=self.number_data
    element = self.element
    label_name = self.label_name
    
    if not os.path.exists(label_dir):
        os.makedirs(label_dir)

    datah = np.zeros((1,2))
    for idx in range(number_data):
        load_bar = (idx/number_data)*(100)
        if load_bar%10 == 0:
            print(str(load_bar) +'% complete')
                
        # write xyz file and convert it to mol file
        at, props = qm9data.get_properties(idx)
        number_atoms = len(props['_atomic_numbers'])
        
        name_xyz = utils.write_xyz_from_db(props,idx,label_dir)
        name_mol = utils.xyz_to_mol(idx,label_dir)
        

        if 1 in props['_atomic_numbers']:        
        #find and store all locations of nitrogens for this molecule
            positions = utils.store_positions(idx,name_mol,props,element)
        
        
        
            # build bond matrix for each of these atoms
            connected_positions = utils.connection_matrix(positions,name_mol,number_atoms)
            
            
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

                nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
                total2,countC2,countH2,countO2,countN2 = utils.count_nn(nn,name_xyz)
                
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
                            
                            array_of_element = utils.np.array([positions[0]])
                            connections = utils.connection_matrix(array_of_element,name_mol,number_atoms)
                            positions2 = connections[0]
                            neighbor2 = min(positions2)

                            index = positions2.index(neighbor2)
                            positions2.pop(index)

                            nn2 = utils.neighboring_connections(name_mol,number_atoms,neighbor2)
                            total3,countC3,countH3,countO3,countN3 = utils.count_nn(nn2,name_xyz)

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
                            connections = utils.connection_matrix(array_of_element,name_mol,number_atoms)
                            positions2 = connections[0]
                            neighbor2 = min(positions2)

                            index = positions2.index(neighbor2)
                            positions2.pop(index)

                            nn2 = utils.neighboring_connections(name_mol,number_atoms,neighbor2)
                            total3,countC3,countH3,countO3,countN3 = utils.count_nn(nn2,name_xyz)

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
                
    savetxt(label_dir+label_name,datah,delimiter=',')
    platform = utils.get_platform()

    if platform == 'Linux':
        os.system('rm ' + label_dir + '*mol')
        os.system('rm ' + label_dir + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_dir = re.sub('../../','/home/mobaxterm/Desktop/DLChem/',label_dir)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_dir + '*xyz"')
