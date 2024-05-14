
'''
Main function that labels dataset according to atom-type (element). Supports H,C,N,O, and functional groups from QM9 
Goes to each atom of a type and determines what functional group is around (through bond connectivity)
'''
# -*- coding: utf-8 -*-

import fileinput
import os
import numpy as np
from schnetpack.datasets import QM9
from numpy import savetxt
import subprocess
import re
from label.manuallabel.utils import utils


class label():
    def __init__(self,qm9data,label_filepath,start,end,element):
        self.label_filepath = label_filepath
        self.start = start
        self.end = end
        self.qm9data = qm9data
        element = ' ' + element + '  '
        self.element = element
        
        if element == ' O  ':
            label_o(self)
        if element == ' N  ':
            raise NotImplementedError
        if element == ' H  ':
            raise NotImplementedError
        if element == ' all  ':
            raise NotImplementedError
        if element == ' C  ':
            raise NotImplementedError
        
def label_o(self):

    allowed_o_labels = {0,3,4,5,7,8,9,10,11,12,14,16,17,21,26,27,28,29}
    acid_importance_list = {8,}

    start = self.start  
    end = self.end
    element = self.element
    label_filepath = self.label_filepath
    qm9data = self.qm9data
    
    label_direc = label_filepath.replace('.csv','/')
    if not os.path.exists(label_direc):
        os.makedirs(label_direc)
    

    datao = np.zeros((1,3))
    countzero = 0


    for idx in range(start,end):
        if idx%1000 == 0:
            print(str(idx) +' complete')
            
        # write xyz file and convert it to mol file
        at, props = qm9data.get_properties(idx)
        number_atoms = len(props['_atomic_numbers'])
        
        name_xyz = utils.write_xyz_from_db(props,idx,label_direc)
        name_mol = utils.xyz_to_mol(idx,label_direc)
        
        if 8 in props['_atomic_numbers']:
            molecule_label_list = [8,16,9,5,4,3]

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
                countF = 0 
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
                            if 'F ' in line2:
                                countF = countF+1
                                
                total = countC + countH + countO + countN+countF
                label = 999
                
                if total == 2:
                    #### OUT OF ORDER IN THE LEGEND INKSCAPE PLOT!!!! BE CAREFUL
                    if countH == 2:
                        label = 16711680
                        marker = 13
                        id = 0

                    if countC == 2:               
                        label = 16747520
                        marker = 7
                        id = 1

                    if countC == 1 and countN == 1:
                        label = 16776960
                        marker = 11
                        id = 2 

                    if countH == 1 and countC == 1:
                        
                        positions = connected_positions[j]
                        neighbor = min(positions)
                        index = positions.index(neighbor)
                        positions.pop(index)
                        
                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
                        total2,countC2,countH2,countO2,countN2,countF2 = utils.count_nn(nn,name_xyz)
                        

                        if total2 == 4: 
                            if countC2 == 3:
                                label = 6737322
                                marker = 1
                                id = 3
                            if countC2 == 2 and countH2 == 1:
                                label = 3978097
                                marker = 15
                                id = 4
                            if countC2 == 1 and countH2 == 2:
                                label = 10025880
                                marker = 9
                                id = 5
                            #### OUT OF ORDER IN THE LEGEND INKSCAPE PLOT!!!! BE CAREFUL
                            if countH2 == 3:
                                label = 10025880
                                marker = 5
                                id = 6
                            #NEWWWWW fooor PKA DATASET 27 & 28!! 
                            if countO2 == 2 and countC2 ==2:
                                label = 13468991
                                marker = 9
                                id = 27
                            if countN2 == 1 and countH2 == 2:
                                label = 5597999
                                marker = 7
                                id = 28
                      
                        if total2 == 3:

                            if countC2 == 2:
                                label = 11393254
                                marker = 2
                                id = 7
                            if countC2 == 1 and countO2 == 2:
                                label = 6266528
                                marker = 5
                                id = 8 
                            if countC2 == 1 and countN2 == 1:
                                label = 8900346
                                marker = 13
                                id = 9 
                            if countO2 == 2 and countN2 == 1:
                                label = 4286945
                                marker = 13
                                id = 10
                            if countN2 == 2:
                                label = 6591981
                                marker = 3
                                id = 11
                            #label 26 is NEW FOR PKA DATA
                            if countC2 == 1 and countH2 == 1:
                                label = 4915330
                                marker = 3
                                id = 26
                            #label 29 is NEW FOR PKA DATA
                            if countO2 == 2 and countH2 == 1:
                                label = 9127187
                                marker = 5
                                id = 29

                         
                    if countH == 1 and countN == 1:
                        label = 15787660
                        marker = 7
                        id = 12
                    if countN == 2:
                        label = 16766720
                        marker = 5
                        id = 13
                if total == 1:
                    if countC == 1:
                        positions = connected_positions[j]
                        neighbor = min(positions)
                        index = positions.index(neighbor)
                        positions.pop(index)
                        
                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
                        total2,countC2,countH2,countO2,countN2, countF2 = utils.count_nn(nn,name_xyz)
                        
                        if total2 == 3: 
                            if countH2 == 2:
                                label = 11674146
                                marker = 15
                                id = 14
                            if countC2 == 2:
                                label = 16758465
                                marker = 11
                                id = 15
                            if countC2 == 1 and countH2 == 1:
                                label = 16770273
                                marker = 1
                                id = 16
                            if countC2 == 1 and countN2 == 1:
                                label = 16716947
                                marker = 3
                                id = 17
                            if countC2 == 1 and countO2 == 2:
                                label = 16738740
                                marker = 15
                                id = 18
                            if countN2 == 1 and countH2 == 1:
                                label = 16716947
                                marker = 3
                                id = 19
                            if countN2 == 1 and countO2 == 2:
                                label = 14524637
                                marker = 5
                                id = 20
                            #label = 25 IS NEW
                            if countN2 == 2 and countO2 == 1:
                                label = 12433259
                                marker = 13
                                id = 25
                            if countN2 == 2:
                                label = 8388736
                                marker = 5
                                id = 21
                            if countO2 == 3:
                                label = 14381203
                                marker = 9
                                id = 22
                            if countO2 == 2 and countH2 == 1:
                                label = 16738740
                                marker = 7
                                id = 23
                    if countN == 1:
                        label = 10824234 
                        marker = 7
                        id = 24
                            
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)

                if id in allowed_o_labels:
                    molecule_label_list.append([label,marker,id,idx])
        print('#labels',len(molecule_label_list))
        if len(molecule_label_list) > 1:
            print(molecule_label_list)
            breakpoint

#                row = np.array([[label,marker,idx]])
#                datao = np.vstack((datao,row))

    datao = np.delete(datao,0,0)                
    savetxt(label_filepath,datao,delimiter=',')
