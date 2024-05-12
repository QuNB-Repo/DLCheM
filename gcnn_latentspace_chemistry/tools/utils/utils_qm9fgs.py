import numpy as np
from numpy import savetxt, genfromtxt
from schnetpack.datasets import QM9
import random

def countfgs(label_filepath, number_fgs,label_column,start_atomidx,end_atomidx,where_is_lonely_class):
    label_data = genfromtxt(label_filepath,delimiter=',',encoding='utf-8-sig',skip_header=1)
    numb_labels = number_fgs

    count_labels = np.zeros(numb_labels)
    for each_atom in range(start_atomidx,end_atomidx):
        label = int(label_data[each_atom][label_column])
        count_labels[label] = count_labels[label] + 1
    
    #if you want to find out which atom_idx 
    if where_is_lonely_class == True:
        for class_ in range(number_fgs):
            if count_labels[class_] == 0:
                print('this class is empty: ', class_)
            if count_labels[class_] == 1:
                atom_idx = np.where(label_data[:,label_column] == class_)
                print('this class has one member: ', class_ , ' and it is located: ', atom_idx)

    print(count_labels)

def extractfg(qm9_filepath,label_filepath,subset_filepath,label_idx):
    qm9data = QM9(qm9_filepath,download=False,remove_uncharacterized=True)
    labeldata = genfromtxt(label_filepath,delimiter=',',encoding='utf-8-sig')
    numb_atoms = len(labeldata)


    #save a list of found fg molec idx to ensure no duplicates are included in the list
    found_mol_list = []
    count_number_molecules = 0

    new_xyz = open(subset_filepath,mode='w')
    for each_atom in range(numb_atoms):

        #ethers
        if labeldata[each_atom][0] == label_idx:

            #find the molecule that atom comes from
            idx = int(labeldata[each_atom][1])

            if idx not in found_mol_list:

                #REMEMBER YOU HAVE TO REMOVE (OR NOT PERMIT) DUPLICATES DUE TO MULTIPLE FGS ON SAME MOL AND COUNT NUMBER OF TOTAL FOUND

                found_mol_list.append(idx)
                count_number_molecules = count_number_molecules + 1
                #retrieve the molecule from qm9data
                at, props = qm9data.get_properties(idx)

                numb_atoms_in_molec = len(props['_positions'])

                new_xyz.write(str(numb_atoms_in_molec) + '\n' + '0.0000' +'\n')
                for each_atom_in_molec in range(numb_atoms_in_molec):
                    if props['_atomic_numbers'][each_atom_in_molec] == 1:
                        new_xyz.write('H ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(props['_positions'][each_atom_in_molec][1]) + ' ' + str(props['_positions'][each_atom_in_molec][2])+'\n')
                    if props['_atomic_numbers'][each_atom_in_molec] == 6:
                        new_xyz.write('C ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(props['_positions'][each_atom_in_molec][1]) + ' ' + str(props['_positions'][each_atom_in_molec][2])+'\n')
                    if props['_atomic_numbers'][each_atom_in_molec] == 7:
                        new_xyz.write('N ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(props['_positions'][each_atom_in_molec][1]) + ' ' + str(props['_positions'][each_atom_in_molec][2])+'\n')
                    if props['_atomic_numbers'][each_atom_in_molec] == 8:
                        new_xyz.write('O ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(props['_positions'][each_atom_in_molec][1]) + ' ' + str(props['_positions'][each_atom_in_molec][2])+'\n')
                    if props['_atomic_numbers'][each_atom_in_molec] == 9:
                        new_xyz.write('F ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(props['_positions'][each_atom_in_molec][1]) + ' ' + str(props['_positions'][each_atom_in_molec][2])+'\n')

    new_xyz.close()
    print(count_number_molecules)



def extractrange(qm9_filepath, subset_filepath, mol_range):
    qm9data = QM9(qm9_filepath, download=False, remove_uncharacterized=True)

    chosen_indices = []  # Store indices of molecules that meet the condition

    full_xyz_file = open(subset_filepath,mode='w')

    for idx in range(mol_range[0],mol_range[1]):
        at, props = qm9data.get_properties(idx)

        if 6 in props['_atomic_numbers']:
            chosen_indices.append(idx)

    count_molecule = 0
    for idx in chosen_indices:

        #Brooke's dataset misses a few molecules! 
        if idx != 45 and idx != 96:
            new_fileend = '%s.xyz' % (str(count_molecule))
            new_xyz = open(subset_filepath.replace('.xyz', new_fileend), mode='w')
            at, props = qm9data.get_properties(idx)

            props['_positions'] = props['_positions'].detach().numpy()

            numb_atoms_in_molec = len(props['_positions'])

            new_xyz.write(str(numb_atoms_in_molec) + '\n' + str(idx) + '\n')
            full_xyz_file.write(str(numb_atoms_in_molec) + '\n' + str(idx) + '\n')

            for each_atom_in_molec in range(numb_atoms_in_molec):
                if props['_atomic_numbers'][each_atom_in_molec] == 1:
                    new_xyz.write(
                        'H ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                    full_xyz_file.write(
                        'H ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                if props['_atomic_numbers'][each_atom_in_molec] == 6:
                    new_xyz.write(
                        'C ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                    full_xyz_file.write(
                        'C ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                if props['_atomic_numbers'][each_atom_in_molec] == 7:
                    new_xyz.write(
                        'N ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                    full_xyz_file.write(
                        'N ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                if props['_atomic_numbers'][each_atom_in_molec] == 8:
                    new_xyz.write(
                        'O ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                    full_xyz_file.write(
                        'O ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                if props['_atomic_numbers'][each_atom_in_molec] == 9:
                    new_xyz.write(
                        'F ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
                    full_xyz_file.write(
                        'F ' + str(props['_positions'][each_atom_in_molec][0]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][1]) + ' ' + str(
                            props['_positions'][each_atom_in_molec][2]) + '\n')
            count_molecule = count_molecule + 1 
            new_xyz.close()
    full_xyz_file.close()

    