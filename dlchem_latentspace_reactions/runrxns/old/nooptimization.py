'''
YOU MAY GET ERRORS IF YOU START IN MIDDLE OF QM9 DATASET
'''

import numpy as np
from numpy import genfromtxt

import fileinput

import os

from fgtransform.utils import utils
from input import *

from schnetpack.datasets import QM9

from rdkit import Chem
from rdkit.Chem import AllChem


#load the molecular dataset that will be used (QM9)
#for the transformation operation
qm9data = QM9(dataset_filepath,download=False,
              remove_uncharacterized=True)

# labelH will be used to identify which
# hydrogen label (aldehyde, alcohol, ... etc)
# to remove in the transformation
labelH = genfromtxt(targetlabelH_filepath,delimiter=' ',encoding='utf-8-sig')


#define the functional group label of H
#that will be used to remove that type of
#H (aldehyde-H, alcohol-H, ...etc)
labelHid = H_id

#define a variable that save the index of H's in label dataset in total, 
#and how many of them are target H found 
n_h_total = 0
n_target = 0

#open dataset in and dataset out
#dataset_in will save all the starting target molecules
#dataset_out will save all the transformed molecules
dataset_in = open(output_filepath1,mode='w')
dataset_out = open(output_filepath2,mode='w')

for idx in range(0,n_molecules):
    #get the molecular signature and properties from 
    #loaded qm9data, define number of atoms
    #from the properties
    at, props = qm9data.get_properties(idx)

    n_atoms = len(props['_atomic_numbers'])

    n_h_molecule = 0 
    for i in range(n_atoms):
        if props['_atomic_numbers'][i] == 1:
            n_h_molecule = n_h_molecule + 1
    
    n_h_total_save = n_h_total
    n_h_total = n_h_total + n_h_molecule 
    
    n_h_fgroup = 0
    for i in range(n_h_total_save,n_h_total):
        if labelH[i][0] == H_id:
            n_h_fgroup = n_h_fgroup + 1


    #only perform operation on molecules with one of this fgroup
    if n_h_fgroup == 1:
        if idx != 80 and idx != 95 and idx != 97:
            print(idx)

            for target in range(n_h_total_save,n_h_total):

                if labelH[target][0] == labelHid: 

                    print(idx)
                    #we found a target! raise the count
                    n_target = n_target + 1

                    #find the hydrogen number of the target within the molecule
                    #take the target number and subtract by the save count of 
                    #total hydrogens just before this molecule
                    #this will be used to identify the hydrogen number in the mol file
                    #to perform the transformation on
                    h_index = target - n_h_total_save

                    #begin writing the xyz file of target molecule
                    fg1tofg2 = ''
                    #xyz file starts with number of atoms at the top
                    fg1tofg2 = fg1tofg2+str(n_atoms) + '\n'
                    #then a name of molecule which we put as 0.0000, no consequence
                    fg1tofg2 = fg1tofg2 + '0.0000 \n'

                    #loop through number of atoms and add their xyz information to string
                    for k in range(n_atoms):
                        if props['_atomic_numbers'][k] == 1:
                            fg1tofg2 = fg1tofg2 + 'H' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]) + '\n'
                        if props['_atomic_numbers'][k] == 6:
                            fg1tofg2 = fg1tofg2 + 'C' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]) + '\n'
                        if props['_atomic_numbers'][k] == 7:
                            fg1tofg2 = fg1tofg2 + 'N' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]) + '\n'
                        if props['_atomic_numbers'][k] == 8:
                            fg1tofg2 = fg1tofg2 + 'O' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]) + '\n'
                        if props['_atomic_numbers'][k] == 9:
                            fg1tofg2 = fg1tofg2 + 'F' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]) + '\n'

                    #open and write xyz file so that it can be converted to mol as RDKit doesn't accept xyz as input
                    #only mol (the program does simple RDKit force field optimization)
                    #also write xyz in the a total target input molecule dataset, dataset_in 
                    init_xyz_filename = 'fg1tofg2-%s-%s.xyz' %(idx,n_target)
                    init_xyz_file = open(init_xyz_filename, mode='w')

                    init_xyz_file.write(fg1tofg2)
                    dataset_in.write(fg1tofg2)
                    init_xyz_file.close()

                    #run obabel to convert xyz to mol
                    init_mol_filename = init_xyz_filename.replace('xyz','mol')
                    os.system('obabel '+init_xyz_filename+' -O '+init_mol_filename)

                    #edit mol file so that C replaces the target H, hydrogen_numb = target - index_h_save 
                    trans_mol_filename = 'fg1tofg2-%s-%strans.mol' %(idx,n_target)
                    trans_mol_file = open(trans_mol_filename,mode='w')
                    count_lines = 0 
                    count_to_h_index = 0 
                    for line in fileinput.FileInput(init_mol_filename, inplace=0):
                        if ' H  ' in line:
                            if count_to_h_index == h_index:
                                line = line.replace('H','C')
                                x = float(line[3:10])
                                y = float(line[13:20])
                                z = float(line[23:30])
                                x = round(1.48*x,4)
                                y = round(1.48*y,4)
                                z = round(1.48*z,4) 
                                if x < 0:
                                    line = line.replace(line[3:10],str(x))
                                if x > 0:
                                    line = line.replace(line[4:10],str(x))
                                if y < 0:
                                    line = line.replace(line[13:20],str(y))
                                if y > 0:
                                    line = line.replace(line[14:20],str(y))
                                if z < 0:
                                    line = line.replace(line[23:30],str(z))
                                if z > 0:
                                    line = line.replace(line[24:30],str(z))

                            count_to_h_index = count_to_h_index + 1
                        trans_mol_file.write(line)

                    trans_mol_file.close()

                    #run RDKit on transformed mol file to add H's
                    molecule = Chem.MolFromMolFile(trans_mol_filename)
                    #add H's to added C
                    molecule_full = Chem.AddHs(molecule)

                    opt = AllChem.MMFFOptimizeMolecule(molecule_full,maxIters=200)

                    molecule_optimized = Chem.MolToMolBlock(molecule_full)

                    final_mol_filename = 'fg1tofg2-%s-%sopt.mol' %(idx,n_target)
                    
                    print('hi',molecule_optimized)

                    opt_file = open(final_mol_filename,mode='w')
                    opt_file.write(molecule_optimized)
                    opt_file.close()

                    #run conversion back to xyz
                    final_xyz_filename = final_mol_filename.replace('mol','xyz')
                    os.system('obabel ' + final_mol_filename + ' -O ' + final_xyz_filename)


                    #write each line of trans xyz file in outputdatset
                    for line in fileinput.FileInput(final_xyz_filename,inplace=0):
                        dataset_out.write(line)
        
print(n_target)
dataset_in.close()
dataset_out.close()

#generate db dataset files from the xyz dataset of initial and transformed
number_inputs = n_target 
available_properties = ['energy']            
utils.generate(output_filepath1,available_properties,number_inputs)
utils.generate(output_filepath2,available_properties,number_inputs)
