# -*- coding: utf-8 -*-
"""
Utility routines for labeller code such as writing xyz data from db file into mol file, convert xyz to mol
with obabel, building bond matrix, storing label_ids, ... etc
"""

import os
import subprocess
import sys
import fileinput
from rdkit import Chem
from rdkit.Chem import AllChem



def xyz2mol(props):

    atomic_nums = props['_atomic_numbers']
    #make xyz block out of molecule for RDKit
    xyz_coords = props['_positions']

    xyz_string = f"{len(atomic_nums)}\n0.000 \n"
    for i in range(len(atomic_nums)):
        xyz_string += f"{atomic_nums[i]} {xyz_coords[i][0]:.6f} {xyz_coords[i][1]:.6f} {xyz_coords[i][2]:.6f}\n"


    #open temporary folder to write xyz into
    file = open('temp.xyz',mode='w')
    file.write(xyz_string)
    file.close()

    os.system('obabel temp.xyz -O temp.mol > /dev/null 2>&1')


#Writes xyz file from a loaded database (db), needs props of db as input, and 
#required for labelling code
def write_xyz_from_db(props,idx,file_path):
    
    name_xyz = file_path+str(idx)+'.xyz'
    # open an empty file
    xyz_file = open(name_xyz,mode='w',encoding='utf-8')

    # copy props['_atomic_numbers'], props['_positions] tensor and change tensor to numpy array 

    #NOTE ON OLDER VERSIONS OF SCHNET YOU HAVE DETACH.NUMPY THE PROPS TENSOR
    atomic_numbers = props['_atomic_numbers'].detach().numpy()
    
    number_atoms = len(atomic_numbers)

    positions = props['_positions'].detach().numpy()
    
    # write xyz file in xyz file format
    xyz_file.write(str(number_atoms)+'\n')
    xyz_file.write('Title'+'\n')
    for i in range(number_atoms):
        if atomic_numbers[i] == 1:
            xyz_file.write('H ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 6:
            xyz_file.write('C ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 7:
            xyz_file.write('N ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 8:
            xyz_file.write('O ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 9:
            xyz_file.write('F ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
    xyz_file.close()
    return name_xyz

#runs obabel conversion tool on the created xyz file to convert to mol file
#required for labelling code
def xyz_to_mol(idx,file_path):
    # define name of temporary xyz file according to idx, define the name of temporary mol file 
    name_xyz = file_path + str(idx) + '.xyz'
    name_mol = file_path + str(idx) + '.mol'
    
    #use obabel to convert xyz to mol
#    output = subprocess.run('obabel ' + name_xyz + ' -O ' + name_mol)

    output = os.system('obabel ' + name_xyz + ' -O ' + name_mol  + ' > /dev/null 2>&1')
    return name_mol



def store_positions(idx,name_mol,props,element):
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

def check_H(neighbor,name_xyz):

    countxyzlines = 0 
    H_present = False
    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'H' in line:
                H_present = True
    return H_present

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

def check_element(neighbor,name_xyz):
    countxyzlines = 0

    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'H ' in line:
                element_present = 'H'
            if 'C ' in line:
                element_present = 'C'
            if 'N ' in line:
                element_present = 'N'
            if 'O ' in line:
                element_present = 'O'
            if 'F ' in line:
                element_present = 'F'

    return element_present

def count_nn(nn,name_xyz):
    countC = 0
    countH = 0
    countO = 0
    countN = 0
    countF = 0
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
                if 'F ' in line2:
                    countF = countF + 1
    total = countC + countH + countO + countN + countF
    
    return total, countC, countH, countO, countN, countF

def create_fglabel(total,countH,countC,countN,countO,countF,branch_id):
    funcgroup_key = ''
    for i in range(countH):
        funcgroup_key = funcgroup_key + 'H'
    for i in range(countC):
        funcgroup_key = funcgroup_key + 'C'
    for i in range(countN):
        funcgroup_key = funcgroup_key + 'N'
    for i in range(countO):
        funcgroup_key = funcgroup_key + 'O'
    for i in range(countF):
        funcgroup_key = funcgroup_key + 'F'    
    funcgroup_key = funcgroup_key + '-'+str(branch_id)+'-'

    return funcgroup_key


