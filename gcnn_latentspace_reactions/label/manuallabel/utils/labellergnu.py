
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
            label_n(self)
        if element == ' H  ':
            label_h(self)
        if element == ' all  ':
            label_all_elements(self)
        if element == ' C  ':
            label_c(self)

def label_all_elements(self):
    start = self.start  
    end = self.end
    label_filepath = self.label_filepath
    qm9data = self.qm9data
    
    label_direc = label_filepath.replace('.csv','')
    if not os.path.exists(label_direc):
        os.makedirs(label_direc)
    
    dataall = np.zeros((1,2))
    for idx in range(start,end):
    
        at, props = qm9data.get_properties(idx)
        
        number_atoms = len(props['_atomic_numbers'])
        atomic_number = props['_atomic_numbers']
        
        for i in range(number_atoms):
            if atomic_number[i] == 1:
                label = 0
                row = [label,idx]
                dataall = np.vstack((dataall,row))
            if atomic_number[i] == 6:
                label = 1
                row = [label,idx]
                dataall = np.vstack((dataall,row))               
            if atomic_number[i] == 7:
                label = 2
                row = [label,idx]
                dataall = np.vstack((dataall,row))
            if atomic_number[i] == 8:
                label = 3
                row = [label,idx]
                dataall = np.vstack((dataall,row))
            if atomic_number[i] == 9:
                label = 4
                row = [label,idx]
                dataall = np.vstack((dataall,row))

    dataall = np.delete(dataall,0,0)
    savetxt(label_filepath,dataall,delimiter=',') 
                
def label_h(self):
    qm9data=self.qm9data
    label_filepath=self.label_filepath
    start = self.start  
    end = self.end
    element = self.element
    
    label_direc = label_filepath.replace('.csv','/')
    if not os.path.exists(label_direc):
        os.makedirs(label_direc)

    datah = np.zeros((1,3))
    countzero = 0
    for idx in range(start,end):
        load_bar = (idx/(end-start))*(100)
        if load_bar%10 == 0:
            print(str(load_bar) +'% complete')
                
        # write xyz file and convert it to mol file
        at, props = qm9data.get_properties(idx)
        number_atoms = len(props['_atomic_numbers'])
        
        name_xyz = utils.write_xyz_from_db(props,idx,label_direc)
        name_mol = utils.xyz_to_mol(idx,label_direc)
        

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
                countF=0
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
                total2,countC2,countH2,countO2,countN2,countF2 = utils.count_nn(nn,name_xyz)
                
                label=999
                if C_nbh == True:
                    if total2 == 4:    
                        if countH2 == 4:
                            #this is methane, no cluster
                            label = 16711680 #red
                            marker = 13     #diamond
                        elif countH2 == 3 and countC2 == 1:
                            #these are all the methyl groups
                            label = 10025880  #palegreen
                            marker = 7       #filled-circle
                        #methyl-O
                        elif countH2 == 3 and countO2 == 1:
                            label = 10145074   #yellowgreen
                            marker = 5         #filled-square
                        #methyl-N
                        elif countH2 == 3 and countN2 == 1:
                            label = 3329330    #limegreen
                            marker = 1         #plus
                        #this is methylene
                        elif countH2 == 2 and countC2 == 2:
                            label = 11529966    #paleturqoise
                            marker = 9         
                        #methylene-N
                        elif countH2==2 and countC2==1 and countN2==1:
                            label = 3978097    #mediumseagreen
                            marker = 1        #
                        elif countH2 == 2 and countC2 == 1 and countO2==1:
                            label = 9109643    #darkmagenta
                            marker = 3         #Asterisk (Star)
                        elif countH2 == 2 and countO2 == 2:
                            label = 4620980    #steelblue
                            marker = 2         #X
                        #methines
                        elif countH2==1 and countC2==3:
                            label=3100495      #darkslategrey
                            marker = 7         #filled-circle
                        elif countH2==1 and countC2==2 and countN2==1:
                            label=35723        #darkcyan
                            marker =11        #down-filled-traingle
                        elif countH2==1 and countC2==2 and countO2==1:
                            label= 4734347      #darkslateblue
                            marker = 2          #X
                        #outer acetal
                        elif countH2==1 and countC2==1 and countO2==2:
                            label=1644912      #midnightblue
                            marker = 3          #Asterisk (star)

                    elif total2 == 3:
                        #trigonal planar carbons
                        #arromatic
                        if countH2==1 and countC2==2:
                            label=9699539       #darkviolet
                            marker= 15          #filled-pentagon
                        elif countH2==1 and countC2==1 and countO2==1:
                            label=8421376       #olive
                            marker=11           #filled-triangle-down
                        elif countH2==1 and countC2==1 and countN2==1:
                            label=4169E1         #royalblue
                            marker = 5           #filled-square
                        elif countH2==1 and countO2==2:
                            label=14315734      
                            marker=15            #pentagon
                        elif countH2==1 and countN2==2:
                            label=9055202
                            marker=1              #plus
                        elif countH2==1 and countN2==1 and countO2==1:
                            label=6266528
                            marker=3              #asterisk/star
                        elif countH2==2 and countO2==1:
                            label=11674146
                            marker=13              #diamon
                            #formaldehyde
                    elif total2 == 2:
                        if countH2==1 and countC2==1:
                            label=32768         #
                            marker=5            #filled-squre
                        elif countH2==1 and countN2==1:
                            label=14204888
                            marker=7            #filled-circle                   
                if N_nbh == True:

                    if total2 == 3:
                        #trigonal planar nitrogens
                        #first one is ammonia, no cluster
                        #ammonia            
                        if countH2==3:
                            label=14423100
                            marker=2
                        #primary amine
                        elif countH2==2 and countC2==1:
                            
                            array_of_element = np.array([positions[0]])
                            connections = utils.connection_matrix(array_of_element,name_mol,number_atoms)
                            positions2 = connections[0]
 
                            check_if_hydrogen = True
                            

                            while check_if_hydrogen == True:
                                neighbor2 = min(positions2)
                                check_if_hydrogen = utils.check_H(neighbor2,name_xyz)
                                index = positions2.index(neighbor2)
                                positions2.pop(index)

                            nn2 = utils.neighboring_connections(name_mol,number_atoms,neighbor2)
                            total3,countC3,countH3,countO3,countN3,countF3 = utils.count_nn(nn2,name_xyz)
                            if total3 == 4:
                                if countC3 == 3:
                                    label=15657130
                                    marker=5
                                if countC3 == 2 and countH3 == 1:
                                    label=13468991
                                    marker=3
                                if countC3 == 1 and countH3 == 2:
                                    label=16768685
                                    marker=5


                            if total3 == 3:
                                if countC3 == 2: 
                                    label=16776960
                                    marker=9
                                if countC3 == 1 and countN3 == 2:
                                    label=16766720
                                    marker=1
                                if countC3 == 1 and countO3 == 1:
                                    label=14329120
                                    marker=3
                                if countO3 == 2:
                                    label=12433259
                                    marker=7
                                if countO3 == 1 and countH3 == 1:
                                    label=16737095
                                    marker=11
                                if countO3 == 1 and countN3 == 2:
                                    label=12092939
                                    marker=15
                                #ANYTHING OUT OF ORDER HAS NOT BEEN ADDED TO COLOR/MARKERS/LEGEND!!!!!!!!!!!!!!!!
                                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                                if countN3 == 2 and countH3 == 1:
#                                    label=6737322
#                                    marker= 1
                                if countN3 == 3:
                                    label=16032864
                                    marker=2

                        #secondary amine
                        elif countH2==1 and countC2==2:
                            label=16753920
                            marker=13
                        elif countH2==1 and countN2==1 and countC2==1:
                            label=16747520
                            marker=1
                        elif countH2==1 and countN2==2:
                            label=13789470
                            marker=15
                    elif total2 == 2:
                        if countH2==1 and countC2==1:
                            label=16770229
                            marker=3
                    elif total2 == 4:
                        if countH2 == 3 and countC2 == 1:
                            label=8388736
                            marker=13
                        if countH2 == 2 and countC2 == 2:
                            label=10824234
                            marker=15
                        if countH2 == 1 and countC2 == 3:
                            label=5597999
                            marker=5
                            
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
                            total3,countC3,countH3,countO3,countN3,countF3 = utils.count_nn(nn2,name_xyz)

### the following was copied, and modified to suit O nb from above, it shows the use in going the extra mile and labelling everything, then you can just copy it!
                            if total3 == 4:
                                if countH3 == 3:
                                    label=8421504
                                    marker=7
                                if countC3 == 3:
                                    label=16767673
                                    marker=13
                                if countC3 == 2 and countH3 == 1:
                                    label=16752762
                                    marker=5
                                if countC3 == 2 and countN3 == 1:
                                    label=16744272
                                    marker=3
                                if countC3 == 2 and countO3 == 2:
                                    label=16416882
                                    marker=2
                                if countC3 == 1 and countH3 == 2:
                                    label=16758465
                                    marker=13
                                    

                            if total3 == 3:
                                if countC3 == 2: 
                                    label=16738740
                                    marker=7
                                if countC3 == 1 and countN3 == 1:
                                    label=16716947
                                    marker=2
                                if countC3 == 1 and countO3 == 2:
                                    label=14381203
                                    marker=1
                                if countO3 == 2 and countN3 == 1:
                                    label=16711935
                                    marker=11
                                if countN3 == 2:
                                    label=15631086
                                    marker=13


                        if countH2 == 2:
                            label=8388608
                            marker=9
                        if countH2==1 and countN2==1:
                            label=16770229
                            marker=3
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)
                        
                row = np.array([[label,marker,idx]])
                datah = np.vstack((datah,row))
    

    datah = np.delete(datah,0,0)
    savetxt(label_filepath,datah,delimiter=',')

def label_c(self):
    qm9data=self.qm9data
    label_filepath=self.label_filepath
    start = self.start  
    end = self.end
    element = self.element
    
    label_direc = label_filepath.replace('.csv','/')
    if not os.path.exists(label_direc):
        os.makedirs(label_direc)

    datac = np.zeros((1,3))
    countzero = 0
    for idx in range(start,end):

        if idx % 1000 == 0:
            print(str(idx) + ' complete')
                
        # write xyz file and convert it to mol file
        at, props = qm9data.get_properties(idx)
        number_atoms = len(props['_atomic_numbers'])
        
        name_xyz = utils.write_xyz_from_db(props,idx,label_direc)
        name_mol = utils.xyz_to_mol(idx,label_direc)
        

        if 6 in props['_atomic_numbers']:        
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
                total = countC + countH + countO + countN + countF

                label=999

                if total == 4:
                    if countC == 4:
                        label=2263842
                        marker=13
                    if countH == 1 and countC == 3:
                        label=10025880
                        marker=7
                    if countH == 2 and countC == 2:
                        label=10145074
                        marker=5
                    if countH == 3 and countC == 1:
                        label=3329330
                        marker=1
                    if countH == 4: 
                        label=16711680
                        marker=9
                    if countH == 2 and countC == 1 and countO == 1:
                        label=11529966
                        marker=1
                    if countH == 2 and countC == 1 and countN == 1:
                        label=9109643
                        marker=3
                    if countH == 3 and countO == 1:
                        label=4620980
                        marker=2
                    if countH == 3 and countN == 1:
                        label=8388736
                        marker=7
                    if countH == 1 and countC == 2 and countO == 1:
                        label=35723
                        marker=11
                    if countH == 1 and countC == 1 and countO == 2:
                        label=6591981
                        marker=2
                    if countH == 1 and countC == 2 and countN == 1:
                        label=1644912
                        marker=3
                    if countC == 3 and countO == 1:
                        label=8900331
                        marker=15
                    if countC == 3 and countN == 1:
                        label=4734347
                        marker=11
                    if countH == 2 and countO == 2:
                        label=4286945
                        marker=5
                    if countF == 4: 
                        label=11674146
                        marker=15
                    if countF == 3 and countC == 1:
                        label=13047173
                        marker=1
                    if countC == 2 and countO == 2:
                        label=2142890
                        marker=3
                if total == 3:
                    if countC == 3:
                        label=9419919
                        marker=13
                    if countH == 1 and countC == 2:
                        label=32768
                        marker=5
                    if countH == 2 and countO == 1:
                        label=14423100
                        marker=7
                    #aldehyde
                    if countH == 1 and countC == 1 and countO == 1:
                        label=2003199
                        marker=13
                    #ketone
                    if countC == 2 and countO == 1:
                        label=3100495
                        marker=5
                    #ester/carboxylicacid
                    if countC == 1 and countO == 2:
                        label=6737322
                        marker=3
                    if countH == 1 and countO == 2:
                        label=11584734
                        marker=5
                    if countH == 1 and countN == 1 and countO ==1:
                        label=16776960
                        marker=9
                    if countC == 1 and countN == 1 and countO ==1:
                        label=16766720
                        marker=1
                    if countC == 1 and countN == 2:
                        label=16711935
                        marker=3
                    if countH == 1 and countC == 1 and countN == 1:
                        label=16738740
                        marker=7
                    if countN == 2 and countO == 1:
                        label=16716947
                        marker=11
                    if countN == 1 and countO == 2:
                        label=14329120
                        marker=15
                    if countN == 3:
                        label=9662683
                        marker=2
                    if countO == 3: 
                        label=16753920
                        marker=13
                    if countC == 2 and countN == 1:
                        label=8087790
                        marker=1
                    if countH == 1 and countN == 2:
                        label=14381203
                        marker=15
                    if countC == 2 and countF == 1:
                        label=16761035
                        marker=3
                    if countC == 1 and countN == 1 and countF == 1:
                        label=12357519
                        marker=13
                    if countN == 2 and countF == 1:
                        label=13458524
                        marker=15
                    
                    
                if total == 2:
                    if countC == 2:
                        label=8421376
                        marker=7
                    if countH == 1 and countC == 1:
                        label=16729344
                        marker=13
                    if countH == 1 and countN == 1:
                        label=14596231
                        marker=5
                    if countC == 1 and countN == 1:
                        label=13468991
                        marker=3
                
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)

                        
                row = np.array([[label,marker,idx]])
                datac = np.vstack((datac,row))

    datac = np.delete(datac,0,0)
    savetxt(label_filepath,datac,delimiter=',')


def label_n(self):
    element = self.element
    start = self.start  
    end = self.end
    label_filepath = self.label_filepath
    qm9data = self.qm9data
    
    label_direc = label_filepath.replace('.csv','/')
    if not os.path.exists(label_direc):
        os.makedirs(label_direc)

    datan = np.zeros((1,3))  
    for idx in range(start,end):
        load_bar = (idx/(end-start))*(100)
        if load_bar%10 == 0:
            print(str(load_bar) +'% complete')
        
        # write xyz file and convert it to mol file
        at, props = qm9data.get_properties(idx)
        number_atoms = len(props['_atomic_numbers'])
        
        name_xyz = utils.write_xyz_from_db(props,idx,label_direc)
        name_mol = utils.xyz_to_mol(idx,label_direc)
        
        
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
                
                label = 999

                if total == 3:
                    if countH == 3:
                        #ammonia NH3
                        label = 16711680
                        marker = 13
                    if countH == 2 and countC == 1:
                        #primary amine H2N-C
                        check_if_hydrogen = True
                        positions = connected_positions[j]
                        while check_if_hydrogen == True:
                            neighbor = min(positions)
                            check_if_hydrogen = utils.check_H(neighbor,name_xyz)
                            index = positions.index(neighbor)
                            positions.pop(index)
                        
                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)


                        total2,countC2,countH2,countO2,countN2, countF2 = utils.count_nn(nn,name_xyz)
                        
                        if total2 == 4: 
                            if countC2 == 3:
                                label=16753920
                                marker=7
                            if countH2 == 1 and countC2 == 2:
                                label=16776960
                                marker=11
                            if countH2 == 2 and countC2 == 1:
                                label=6737322
                                marker=1
                                                            
                        if total2 == 3:
                            if countC2 == 2:
                                label=3978097
                                marker=15
                            if countC2 == 1 and countO2 == 1:
                                label=10025880
                                marker=9
                            if countO2 == 2:
                                label=10145074
                                marker=5
                            if countO2 == 1 and countN2 == 2:
                                label=11393254
                                marker=15
                            if countN2 == 3: 
                                label=6266528
                                marker=3
                            if countN2 == 2 and countC2 == 1:
                                label=8900346
                                marker=13
                            if countN2 == 2 and countH2 == 1:
                                label=4286945
                                marker=13
                            if countO2 == 1 and countH2 ==1:
                                label=6591981
                                marker=3
                
                    if countH == 1 and countC == 2:
                        label=15787660
                        marker=7
                    if countH == 1 and countN == 2:
                        label=16766720
                        marker=5
                    if countH ==1 and countN == 1 and countC ==1:
                        label=13468991
                        marker=15
                    if countC == 3:
                        label=16758465
                        marker=11
                    if countC == 2 and countN == 1:
                        label=16770273
                        marker=1
                    if countC == 1 and countN ==2:
                        label=16716947
                        marker=3
                if total == 2:
                    if countC == 2:
                        label=16738740
                        marker=15
                    if countC == 1 and countH == 1:
                        label=13047173
                        marker=1
                    if countC == 1 and countO == 1:
                        label=14524637
                        marker=13
                    if countC == 1 and countN == 1:
                        label=8388736
                        marker=5
                    if countO == 1 and countN == 1:
                        label=14381203
                        marker=9
                    if countN == 2:
                        label=12211667
                        marker=7
                if total == 1:
                    if countC == 1:
                        label=13808780
                        marker=5
                if total == 4:
                    if countH == 3 and countC == 1:
                        label=13458524
                        marker=7
                    if countH == 2 and countC == 2:
                        label=16752762
                        marker=9
                    if countC == 3 and countH == 1:
                        label=15761536
                        marker=2
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)

                        
                row = np.array([[label,marker,idx]])
                datan = np.vstack((datan,row))

    datan = np.delete(datan,0,0)
    savetxt(label_filepath,datan,delimiter=',')
        
def label_o(self):

    allowed_o_labels = {0,3,4,5,7,8,9,10,11,12,14,16,17,21,26,27,28,29}

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
                                label = 14423100
                                marker = 5
                                id = 6
                            #NEWWWWW fooor PKA DATASET
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
                                id = 17
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
                                id = 18
                    if countN == 1:
                        label = 10824234 
                        marker = 7
                        id = 24
                            
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)

                row = np.array([[label,marker,idx]])
                datao = np.vstack((datao,row))

    datao = np.delete(datao,0,0)                
    savetxt(label_filepath,datao,delimiter=',')
