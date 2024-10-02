'''
PLEASE NOTE YOU CANNOT USE MIDDLE OF DATASET
MUST BEGIN AT THE START OF DATASET

inputs: qm9data, embedding data, >

outputs: alcohol embedding, const>

This program investigates if ther>
Can the embedding vector of CH3OC>
'''

#from distutils.command.check import HAS_DOCUTILS
#from unittest.util import _count_diff_hashable

#from torch import hann_window
import numpy as np
from numpy import genfromtxt
from schnetpack.datasets import QM9
import fileinput
import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from input import *
from fgtransform.utils import utils

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
index_h = 0
count_target = 0

#open dataset in and dataset out
#dataset_in will save all the starting target molecules
#dataset_out will save all the transformed molecules
dataset_in = open(output_filepath1,mode='w')
dataset_out = open(output_filepath2,mode='w')

#run through a certain number of QM9 molecules in dataset
#CANNOT START IN THE MIDDLE
for idx in range(0,n_molecules):


    #get the molecular signature and properties from 
    #loaded qm9data, define number of atoms
     #from the properties
    at, props = qm9data.get_properties(idx)


    props['_atomic_numbers'] = props['_atomic_numbers'].detach().numpy()
    props['_positions'] = props['_positions'].detach().numpy()

    numb_atoms = len(props['_atomic_numbers'])

    #check if this molecule contains one and ONLY ONE oxygen
    number_element = 0
    numb_h = 0 
    for i in range(numb_atoms):
        if props['_atomic_numbers'][i] == atom_num:
            number_element = number_element + 1

        #begin a couter of the number of hydrogens
        #in this molecule alone
        if props['_atomic_numbers'][i] == 1:
            numb_h = numb_h + 1

    #save the initial index of hydrogen that you are on in the label file
    index_h_save = index_h 

    #save the new hydrogen index, since it is index it should be subtracted by 1
    #we should save index as later we will use the hydrogen index as a
    #counting minimum towards the target hydrogen in the molecule,
    #use the total count_h
    index_h = index_h + numb_h 


    #ONLY include those molecules where there is only one target element (this cannot be done with C's and H's)
    if number_element == 1: 
            

    #    aldskets
    #    if idx != 59 and idx != 189 and idx != 190 and idx != 191 and idx != 201 and idx != 202 and idx != 205 and idx != 239 and idx != 255 and idx != 257 and idx != 334 and idx != 744 and idx != 748 and idx != 749 and idx != 756:

    #    1alcseths
    #    print(idx)
    #    if idx != 76 and idx != 262 and idx != 369 and idx != 540 and idx != 1136 and idx != 1138 and idx != 1515 and idx != 4687 and idx != 4688 and idx != 4749 and idx != 4751 and idx != 4806:

    #    2alcseths
    #    if idx != 371 and idx != 1579 and idx != 1719 and idx != 3031:
    #

    #    3alcseths
    #    if idx != 1139 and idx != 1721:   
    #    1amns2amns
    #    if idx != 370 and idx != 815 and idx != 1137 and idx != 4691 and idx != 4692:
    #     
 #        2amns?
 #       if idx != 861 and idx != 957 and idx != 4350 and idx != 4368 and idx != 4380 and idx != 4390 and idx != 4422 and idx != 4945:
 #        1amnseths
#        if idx != 4691: 
# 
#       HON
#        if idx != 3977 and idx != 4028:

#new code
#ohch2ch2c-ch3och2ch2c
        if idx != 1515 and idx != 5797 and idx != 9122:
            print(number_element)
            for target in range(index_h_save,index_h):

                if labelH[target][0] == labelHid:
                    
                    print(idx)
                    #we found a target! raise the count
                    count_target = count_target + 1

                    #find the hydrogen number of the target within the molecule
                    #take the target number and subtract by the save count of 
                    #total hydrogens just before this molecule
                    #this will be used to identify the hydrogen number in the mol file
                    #to perform the transformation on
                    hydrogen_numb = target - index_h_save

                    #begin writing the xyz file of target molecule
                    fg1tofg2 = ''
                    #xyz file starts with number of atoms at the top
                    fg1tofg2 = fg1tofg2+str(numb_atoms) + '\n'
                    #then a name of molecule which we put as 0.0000, no consequence
                    fg1tofg2 = fg1tofg2 + '0.0000 \n'

                    #loop through number of atoms and add their xyz information to string
                    for k in range(numb_atoms):
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
                    init_xyz_filename = 'fg1tofg2-%s-%s.xyz' %(idx,count_target)
                    init_xyz_file = open(init_xyz_filename, mode='w')

                    init_xyz_file.write(fg1tofg2)
                    dataset_in.write(fg1tofg2)
                    init_xyz_file.close()

                    #run obabel to convert xyz to mol
                    init_mol_filename = init_xyz_filename.replace('xyz','mol')
                    os.system('obabel '+init_xyz_filename+' -O '+init_mol_filename)

                    #edit mol file so that C replaces the target H, hydrogen_numb = target - index_h_save 
                    trans_mol_filename = 'fg1tofg2-%s-%strans.mol' %(idx,count_target)
                    trans_mol_file = open(trans_mol_filename,mode='w')
                    count_lines = 0 
                    count_to_hydrogen_numb = 0 
                    for line in fileinput.FileInput(init_mol_filename, inplace=0):
                        if ' H  ' in line:
                            if count_to_hydrogen_numb == hydrogen_numb:
                                line = line.replace('H','C')

                            count_to_hydrogen_numb = count_to_hydrogen_numb + 1
                        trans_mol_file.write(line)

                    trans_mol_file.close()

                    #run RDKit on transformed mol file to optimize it with MMFF force field
                    molecule = Chem.MolFromMolFile(trans_mol_filename)
                    #add H's to added C
                    molecule_full = Chem.AddHs(molecule)
                    #optimize full molecule with added Hs
                    opt = AllChem.MMFFOptimizeMolecule(molecule_full)

                    molecule_optimized = Chem.MolToMolBlock(molecule_full)

                    final_mol_filename = 'fg1tofg2-%s-%sopt.mol' %(idx,count_target)
                    
                    opt_file = open(final_mol_filename,mode='w')
                    opt_file.write(molecule_optimized)
                    opt_file.close()

                    #run conversion back to xyz
                    final_xyz_filename = final_mol_filename.replace('mol','xyz')
                    os.system('obabel ' + final_mol_filename + ' -O ' + final_xyz_filename)

                    for line in fileinput.FileInput(final_xyz_filename,inplace=0):
                        dataset_out.write(line)
        
print(count_target)
dataset_in.close()
dataset_out.close()

#generate db dataset files from the xyz dataset of initial and transformed
number_inputs = count_target 
available_properties = ['energy']            
utils.generate(output_filepath1,available_properties,number_inputs)
utils.generate(output_filepath2,available_properties,number_inputs)
