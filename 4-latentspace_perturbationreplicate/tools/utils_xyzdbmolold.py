import numpy as np
from numpy import savetxt, genfromtxt
from ase.io import read
from schnetpack import AtomsData
from schnetpack.datasets import QM9
import fileinput
import os

def atomic_number_to_element_symbol(atomic_id):
    element_symbols = {
            1: "H",
            6: "C",
            7: "N",
            8: "O",
            9: "F"
        }
    
#    if reverse == False: 
#        return element_symbols.get(atomic_id, None)

#    if reverse == True:
#        reversed_dict = {}
#        for key, value in element_symbols.items():
#            reversed_dict[value] = key

    return element_symbols.get(atomic_id,None)

def atomic_xyz_string(atomic_numbers, atomic_positions):
    if len(atomic_numbers) != len(atomic_positions):
        raise ValueError("The number of atomic numbers and atomic positions should be the same.")

    element_symbols = [atomic_number_to_element_symbol(num) for num in atomic_numbers]
    num_atoms = len(atomic_numbers)
    xyz_string = f"{num_atoms}\n0.00000\n"

    for element_symbol, (x, y, z) in zip(element_symbols, atomic_positions):
        xyz_string += f"{element_symbol} {x:.5f} {y:.5f} {z:.5f}\n"

    return xyz_string


def bunch_xyz(xyz_folder_dir,big_xyzfilename):

    big_xyzfile = open(big_xyzfilename,mode='w')
    for filename in os.listdir(xyz_folder_dir):
        if os.path.isfile(os.path.join(xyz_folder_dir,filename)):
            if filename.endswith('.xyz'):

                xyz_file = open(xyz_folder_dir+filename,mode='r')

                for index, line in enumerate(xyz_file.read().split('\n')):
                    
                    if index == 1:
                        line = line.replace(line,'0.0000')
                    
                    if line != '':
                        big_xyzfile.write(line+'\n')
    
                xyz_file.close()

    big_xyzfile.close()                



def gendb(datasets_filepath,available_properties,index):

    index = str(index)
    atoms = read(datasets_filepath, index=':'+index)
    
    property_list = []
    
    for at in atoms:
        energy = np.array([float(list(at.info.keys())[0])], dtype=np.float32)
        property_list.append(
            {'energy': energy}
        )
        
    print('Properties:', property_list)
    
    new_dataset_filepath = datasets_filepath.replace('.xyz','.db')
    new_dataset = AtomsData(new_dataset_filepath, available_properties=[available_properties])
    new_dataset.add_systems(atoms,property_list)

def countxyz(filepath):
    count = 0
    for line in fileinput.FileInput(filepath,inplace=0):
        if '0.0000' in line:
            count = count + 1
    print(count)

def usedb(db_filepath,idx,qm9=True,available_properties=['NA']):
    if qm9 == True:
        data = QM9(db_filepath,download=False,remove_uncharacterized=True)
        at, props = data.get_properties(idx)

        props['_positions'] = props['_positions'].detach().numpy()
        props['_atomic_numbers'] = props['_atomic_numbers'].detach().numpy()
        
        for k in range(len(props['_atomic_numbers'])):
            if props['_atomic_numbers'][k] == 1:
                print('H' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 6:
                print('C' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 7:
                print('N' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 8:
                print('O' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 9:
                print('F' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))

    else:
        data = AtomsData(db_filepath,available_properties=available_properties)
    
        at, props = data.get_properties(idx)

        props['_positions'] = props['_positions'].detach().numpy()
        props['_atomic_numbers'] = props['_atomic_numbers'].detach().numpy()
        

        for k in range(len(props['_atomic_numbers'])):
            if props['_atomic_numbers'][k] == 1:
                print('H' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 6:
                print('C' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 7:
                print('N' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 8:
                print('O' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))
            if props['_atomic_numbers'][k] == 9:
                print( 'F' + ' ' + str(props['_positions'][k][0]) + ' ' + str(props['_positions'][k][1]) + ' ' + str(props['_positions'][k][2]))


def convertxyz_to_mol(input_xyz_filename,output_mol_filename):

    out = os.system('obabel ' + input_xyz_filename +' -O ' + output_mol_filename + '')

def convertxyz_to_molfull(atomicnumbers,atomicpositions,temp_filename):
    molxyz = atomic_xyz_string(atomicnumbers,atomicpositions)

    tempxyzfilepath = temp_filename+".xyz"
    with open(tempxyzfilepath, 'w') as file:
        file.write(molxyz)

    tempmolfilepath = temp_filename+".mol"
    convertxyz_to_mol(tempxyzfilepath,tempmolfilepath)

    return tempmolfilepath

def countatomsinmol(mol):

    numberatoms = 0 
    for line_index, each_line in enumerate(mol):
        if line_index > 3:
            if 'H' or ' C ' or ' O ' or ' N ' or ' F ' in each_line:
                numberatoms = numberatoms+1 
    
    numberbonds = 0
    for line_index, each_line in enumerate(mol):
        if line_index > 3 + numberatoms:
            if 'END' not in each_line:
                numberbonds = numberbonds + 1


    return numberatoms,numberbonds