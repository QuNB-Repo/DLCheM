from input2 import *

from schnetpack.datasets import QM9

import numpy as np
from numpy import genfromtxt, savetxt

import os
import fileinput

from fgtransform.utils import utils

#load dataset of molecules, db file
qm9data = QM9(dataset_filepath,download=False, remove_uncharacterized=True)

#load dataset labels for H and C
labelH = genfromtxt(targetlabelH_filepath,delimiter=',',encoding='utf-8-sig')

#rename target label id
labelH1id = H1_id
labelH2id = H2_id

#count_h
#count_target
index_h = 0 
count_target = 0 

#initiate init and trans datasets
dataset_in = open(output_filepath1, mode='w')
dataset_out = open(output_filepath2, mode='w')



#for range of all molecules starting from 0
for idx in range(0,n_molecules):



    #get the properties of molecule
    at, props = qm9data.get_properties(idx)

    #get the numb of atoms
    numb_atoms = len(props['_atomic_numbers'])

    #initiate number of target element per molecule (to check there is only one per molecule)
    number_element = 0 
    #initiate number of H's
    numb_h = 0

    for i in range(numb_atoms):
        #count number of target element per molecule
        if props['_atomic_numbers'][i] == atom_num:
            number_element = number_element + 1
        
        #count the number of hydrogens
        if props['_atomic_numbers'][i] == 1:
            numb_h = numb_h + 1
    
    #save previous h index
    index_h_save = index_h

    #save new h index
    index_h = index_h + numb_h


    #only do the following if the number of target element found is 1
    if number_element == 1:
        #target1, target2 switches for checking both are satisfied
        target1 = False
        target2 = False

        #check between the saved H indices for the target H label
        for target in range(index_h_save,index_h):


            #check if target H index in label has the right target H ID 1
            if labelH[target][0] == labelH1id:
                target1 = True


                #count targets found
                count_target = count_target + 1

                #save hydrogen number in the molecule for use in mol file search, target - previous_index
                hydrogen_numb = target - index_h_save

            #check if target H index in label has the right target H ID 2
            if labelH[target][0] == labelH2id:
                target2 = True

                #save hydrogen number in the molecule for use in mol file search, target - previous_index
                hydrogen_numb2 = target - index_h_save
        
        
        if target1 == True and target2 == True:
            #write the xyz of molecule
            fg1tofg2 = ''

            fg1tofg2 = fg1tofg2 + str(numb_atoms) + '\n'

            fg1tofg2 = fg1tofg2 + '0.0000 \n'
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

            #open xyz file
            init_xyz_filename = 'fg1tofg2-%s-%s.xyz' %(idx,count_target)
            init_xyz_file = open(init_xyz_filename, mode='w')
                
            #write molecule on it
            init_xyz_file.write(fg1tofg2)

            #write init dataset
            dataset_in.write(fg1tofg2)

            #close xyz file
            init_xyz_file.close()

            #make the oxidation edits! 
                
            #count lines
            count_lines = 0 
            count_atoms = 0
            count_to_hydrogen_numb = 0

            #for each line in init mol file
            trans_xyz_string = ''
            for line in fileinput.FileInput(init_xyz_filename,inplace=0):
                count_lines = count_lines + 1
                
                if count_lines > 2:
                    count_atoms = count_atoms + 1

                    #if you find an H in line
                    if 'H' in line:
                        #count to the hydrogen number of target
                        if count_to_hydrogen_numb == hydrogen_numb:
                            line = line.replace(line,'')
                            count_to_h1 = count_atoms 

                        if count_to_hydrogen_numb == hydrogen_numb2:
                            line = line.replace(line,'')
                            count_to_h2 = count_atoms

                        count_to_hydrogen_numb = count_to_hydrogen_numb + 1
                    
                trans_xyz_string = trans_xyz_string +  line 

            #final touches, change atom number, -2 , bond number, -2, remove empty lines after line 4,
            #add double bond in the right place

            #open trans mol file name
            trans_xyz_filename = 'fg1tofg2-%s-%strans.xyz' %(idx,count_target)
            trans_xyz_file = open(trans_xyz_filename,mode='w')
            
            count_lines = 0 
            trans=''
            for line in trans_xyz_string.split('\n'):
                count_lines = count_lines + 1

                #change atom  number in xyz file -2
                if count_lines == 1:
                    line = line.replace(line,str(int(line)-2))

                if line != '\n' and count_lines < 1 + numb_atoms:
                    trans = trans + line + '\n'
                    trans_xyz_file.write(line + '\n') 

                if count_lines == 1 + numb_atoms: 
                    trans = trans + line
                    trans_xyz_file.write(line)   


            trans_xyz_file.close()

            dataset_out.write(trans)


dataset_in.close()
dataset_out.close()


number_inputs = count_target
print('number_inputs',number_inputs) 
available_properties = ['energy']            
utils.generate(output_filepath1,available_properties,number_inputs)
utils.generate(output_filepath2,available_properties,number_inputs)






