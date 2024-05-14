
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


    platform = utils.get_platform()
    
    #rtied a million times not to use MOBA... there is no combatibility. MOBA must be used at the end because it needs to be closed
    #has to do with the backslash in windows vs how it is used here
    if platform == 'Linux':
        os.system('rm ' + label_direc + '*mol')
        os.system('rm ' + label_direc + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_direc = re.sub('../../','/home/mobaxterm/Desktop/DLChem/',label_direc)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*xyz"')    
                
        
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
                            marker = 9         #up-filled-triangle
                        #methylene-N
                        elif countH2==2 and countC2==1 and countN2==1:
                            label = 3978097    #mediumseagreen
                            marker = 15        #filled-pentagon
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
                            marker = 11        #down-filled-traingle
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
                            marker=13             #diamon
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
                                if countN3 == 2 and countH3 == 1:
                                    label=6737322
                                    marker= 7
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
                            label=15787660
                            marker=15
                        if countH2 == 2 and countC2 == 2:
                            label=8388736
                            marker=15
                        if countH2 == 1 and countC2 == 3:
                            label=10824234
                            marker=7
                            
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
                                    marker=15
                                if countC3 == 3:
                                    label=16767673
                                    marker=5
                                if countC3 == 2 and countH3 == 1:
                                    label=16752762
                                    marker=3
                                if countC3 == 2 and countN3 == 1:
                                    label=16744272
                                    marker=2
                                if countC3 == 2 and countO3 == 2:
                                    label=16416882
                                    marker=13
                                if countC3 == 1 and countH3 == 2:
                                    label=16758465
                                    marker=7
                                    

                            if total3 == 3:
                                if countC3 == 2: 
                                    label=16738740
                                    marker=2
                                if countC3 == 1 and countN3 == 1:
                                    label=16716947
                                    marker=15
                                if countC3 == 1 and countO3 == 2:
                                    label=14381203
                                    marker=11
                                if countO3 == 2 and countN3 == 1:
                                    label=16711935
                                    marker=13
                                if countN3 == 2:
                                    label=15631086
                                    marker=9


                        if countH2 == 2:
                            label=8388608
                            marker=2
                        if countH2==1 and countN2==1:
                            label=16770229
                            marker=7
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)
                        
                row = np.array([[label,marker,idx]])
                datah = np.vstack((datah,row))
    

    datah = np.delete(datah,0,0)
    savetxt(label_filepath,datah,delimiter=',')
    platform = utils.get_platform()

    if platform == 'Linux':
        os.system('rm ' + label_direc + '*mol')
        os.system('rm ' + label_direc + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        subprocess.run(f'{moba_path} -newtab "rm /home/mobaxterm/Desktop/Amer/UNB/Projects/DLChem/data/labelqm9/H/*mol')
        subprocess.run(f'{moba_path} -newtab "rm /home/mobaxterm/Desktop/Amer/UNB/Projects/DLChem/data/labelqm9/H/*xyz')


def label_c(self):
    qm9data=self.qm9data
    label_filepath=self.label_filepath
    start = self.start  
    end = self.end
    element = self.element
    
    label_direc = label_filepath.replace('.csv','/')
    if not os.path.exists(label_direc):
        os.makedirs(label_direc)

    datac = np.zeros((1,2))
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
                        label = 0
                    if countH == 1 and countC == 3:
                        label = 1
                    if countH == 2 and countC == 2:
                        label = 2
                    if countH == 3 and countC == 1:
                        label = 3
                    if countH == 4: 
                        label = 4
                    if countH == 2 and countC == 1 and countO == 1:
                        label = 5
                    if countH == 2 and countC == 1 and countN == 1:
                        label = 6
                    if countH == 3 and countO == 1:
                        label = 7
                    if countH == 3 and countN == 1:
                        label = 8
                    if countH == 1 and countC == 2 and countO == 1:
                        label = 9
                    if countH == 1 and countC == 1 and countO == 2:
                        label = 10
                    if countH == 1 and countC == 2 and countN == 1:
                        label = 11
                    if countC == 3 and countO == 1:
                        label = 12
                    if countC == 3 and countN == 1:
                        label = 13
                    if countH == 2 and countO == 2:
                        label = 14
                    if countF == 4: 
                        label = 15
                    if countF == 3 and countC == 1:
                        label = 16
                    if countC == 2 and countO == 2:
                        label = 17
                if total == 3:
                    if countC == 3:
                        label = 18
                    if countH == 1 and countC == 2:
                        label = 19
                    if countH == 2 and countO == 1:
                        label = 20
                    #aldehyde
                    if countH == 1 and countC == 1 and countO == 1:
                        label = 21
                    #ketone
                    if countC == 2 and countO == 1:
                        label = 22
                    #ester/carboxylicacid
                    if countC == 1 and countO == 2:
                        label = 23
                    if countH == 1 and countO == 2:
                        label = 24
                    if countH == 1 and countN == 1 and countO ==1:
                        label = 25
                    if countC == 1 and countN == 1 and countO ==1:
                        label = 26
                    if countC == 1 and countN == 2:
                        label = 27
                    if countH == 1 and countC == 1 and countN == 1:
                        label = 28
                    if countN == 2 and countO == 1:
                        label = 29
                    if countN == 1 and countO == 2:
                        label = 30
                    if countN == 3:
                        label = 31
                    if countO == 3: 
                        label = 32
                    if countC == 2 and countN == 1:
                        label = 33
                    if countH == 1 and countN == 2:
                        label = 34
                    if countC == 2 and countF == 1:
                        label = 35
                    if countC == 1 and countN == 1 and countF == 1:
                        label = 36
                    if countN == 2 and countF == 1:
                        label = 37
                    
                    
                if total == 2:
                    if countC == 2:
                        label = 38
                    if countH == 1 and countC == 1:
                        label = 39
                    if countH == 1 and countN == 1:
                        label = 40
                    if countC == 1 and countN == 1:
                        label = 41
                
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)

                        
                row = np.array([[label,idx]])
                datac = np.vstack((datac,row))

    datac = np.delete(datac,0,0)
    savetxt(label_filepath,datac,delimiter=',')

    platform = utils.get_platform()
    
    if platform == 'Linux':
        os.system('rm ' + label_direc + '*mol')
        os.system('rm ' + label_direc + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_direc = re.sub('../../','/home/mobaxterm/Desktop/DLChem/',label_dir)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*xyz"')
        
             


def label_n(self):
    element = self.element
    start = self.start  
    end = self.end
    label_filepath = self.label_filepath
    qm9data = self.qm9data
    
    label_direc = label_filepath.replace('.csv','/')
    if not os.path.exists(label_direc):
        os.makedirs(label_direc)

    datan = np.zeros((1,2))  
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
                total = countC + countH + countO + countN + countF
                
                label = 999

                if total == 3:
                    if countH == 3:
                        #ammonia NH3
                        label = 0
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


                        total2,countC2,countH2,countO2,countN2,countF2 = utils.count_nn(nn,name_xyz)
                        
                        if total2 == 4: 
                            if countC2 == 3:
                                label = 1
                            if countH2 == 1 and countC2 == 2:
                                label = 2
                            if countH2 == 2 and countC2 == 1:
                                label = 3
                                                            
                        if total2 == 3:
                            if countC2 == 2:
                                label = 4
                            if countC2 == 1 and countO2 == 1:
                                label = 5
                            if countO2 == 2:
                                label = 6
                            if countO2 == 1 and countN2 == 2:
                                label = 7
                            if countN2 == 3: 
                                label = 8
                            if countN2 == 2 and countC2 == 1:
                                label = 9
                            if countN2 == 2 and countH2 == 1:
                                label = 10
                            if countO2 == 1 and countH2 ==1:
                                label = 11
                
                    if countH == 1 and countC == 2:
                        label = 12    
                    if countH == 1 and countN == 2:
                        label = 13
                    if countH ==1 and countN == 1 and countC ==1:
                        label = 14
                    if countC == 3:
                        label=15
                    if countC == 2 and countN == 1:
                        label = 16
                    if countC == 1 and countN ==2:
                        label = 17
                if total == 2:
                    if countC == 2:
                        label = 18
                    if countC == 1 and countH == 1:
                        label = 19
                    if countC == 1 and countO == 1:
                        label = 20
                    if countC == 1 and countN == 1:
                        label= 21
                    if countO == 1 and countN == 1:
                        label = 22 
                    if countN == 2:
                        label = 23
                if total == 1:
                    if countC == 1:
                        label = 24
                if total == 4:
                    if countH == 3 and countC == 1:
                        label = 25
                    if countH == 2 and countC == 2:
                        label = 26
                    if countC == 3 and countH == 1:
                        label = 27
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)

                        
                row = np.array([[label,idx]])
                datan = np.vstack((datan,row))

    datan = np.delete(datan,0,0)
    savetxt(label_filepath,datan,delimiter=',')

    platform = utils.get_platform()
    
    if platform == 'Linux':
        os.system('rm ' + label_direc + '*mol')
        os.system('rm ' + label_direc + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_direc = re.sub('../../','/home/mobaxterm/Desktop/DLChem/',label_dir)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*xyz"')
        

        
def label_o(self):
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
                label = 999
                
                if total == 2:
                    if countH == 2:
                        id = 0 
                        label = 16711680
                        marker = 13

                    if countC == 2:
                        id = 1                
                        label = 16747520
                        marker = 7

                    if countC == 1 and countN == 1:
                        label = 16776960
                        marker = 11
                        id = 2

                    if countH == 1 and countC == 1:
                        label = 11403055
                        marker = 15
                        id = 3

                    if countH == 1 and countN == 1:
                        label = 15787660
                        marker = 7
                        id = 4

                    if countN == 2:
                        label = 16766720
                        marker = 5
                        id = 5
                if total == 1:
                    if countC == 1:
                        label = 14381203
                        marker = 1
                        id = 6
                    if countN == 1:
                        label = 10824234
                        marker = 7
                        id = 7
                            
                if label == 999:
                    print('999')
                    error_message = 'this functional group is unknown, molecule_id = %s' %(idx)
                    raise ValueError(error_message)


                row = np.array([[label,marker,idx]])
                datao = np.vstack((datao,row))

    datao = np.delete(datao,0,0)                
    savetxt(label_filepath,datao,delimiter=',')


    platform = utils.get_platform()

    #rtied a million times not to use MOBA... there is no combatibility. MOBA must be used at the end because it needs to be closed
    #has to do with the backslash in windows vs how it is used here
    if platform == 'Linux':
        os.system('rm ' + label_direc + '*mol')
        os.system('rm ' + label_direc + '*xyz')
    if platform == 'Windows':
        moba_path = "C:\Program Files (x86)\Mobatek\MobaXterm\MobaXterm.exe"
        label_direc = re.sub('../../','/home/mobaxterm/Desktop/DLChem/',label_direc)
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*mol"')
        subprocess.run(f'{moba_path} -newtab "rm ' + label_direc + '*xyz"')

