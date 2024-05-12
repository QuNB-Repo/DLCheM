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
    

def generate_dbfromxyz(datasets_filepath,available_properties):
    '''
    generate db file from xyz database file 
    xyz database file will should look like this:

            #atoms
            property1 property2 ...
            element x y z 
            element x y z
            ...
            #atoms
            property1 property2 ...
            element x y z
            element x y z 
            ...
            ...
    
        Args:
            dataset_filepath            - the xyz dataset filepath to be converted to db
            available_properties        - list of strings representing the names of the properties
            number_properties           - number of properties for each molecule

        Process:
            molecules                   - uses the atomic simulation environment (ase) package to read 
                                          the xyz dataset filepath and ensure that it is separated and indexed
                                          by molecular 'frames'
            property_list               - list that will hold all the properties 
            mol                         - each molecule in the xyz ase read file
            each_prop                   - index that runs over all the properties available

        
        Returns:
            a db file path next to the xyz file it was converted from 
            db files can be accessed with SchNet using AtomsData
    '''

    molecules = read(datasets_filepath, index=':')
    
    number_properties  = len(available_properties)

    properties_list = []

    for mol in molecules:

        print(mol.info.keys())

        properties = [np.array([float(list(mol.info.keys())[each_prop])], dtype=np.float32) for each_prop in range(number_properties)]
        print(properties)

        properties_list.append({available_properties[each_prop]: properties[each_prop] for each_prop in range(number_properties)})

    print('Properties:', properties_list)
    
    
    new_dataset_filepath = datasets_filepath.replace('.xyz','.db')
    new_dataset = AtomsData(new_dataset_filepath, available_properties=available_properties)
    new_dataset.add_systems(molecules,properties_list)



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
    
    os.system('obabel ' + input_xyz_filename +' -O ' + output_mol_filename )

def convertxyz_to_molfull(atomicnumbers,atomicpositions):
    molxyz = atomic_xyz_string(atomicnumbers,atomicpositions)

    tempxyzfilepath = "temp.xyz"
    with open(tempxyzfilepath, 'w') as file:
        file.write(molxyz)

    tempmolfilepath = "temp.mol"
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