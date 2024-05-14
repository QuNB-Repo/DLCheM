import os

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

    os.system('obabel temp.xyz -O temp.mol')
    
    
    mol_filename = 'temp.mol'
    return mol_filename


def find_connections(number_atoms,each_atom,mol_file_read):

    connections_to_atom = []

    #go to the bonding part of the molfile, find what connects to each atom
    for line_index, each_line in enumerate(mol_file_read.split('\n')):
        if line_index > number_atoms + 3:
            #little line processing to remove nasty repeats
            each_line = each_line.replace(each_line,each_line[:7])
            if ' '+ str(each_atom)+' ' in each_line and 'RAD' not in each_line:
                if int(each_line[1:4]) == int(each_atom) or int(each_line[0:4]) == int(each_atom):
                    connections_to_atom.append(int(each_line[4:6]))
                else:
                    connections_to_atom.append(int(each_line[1:3]))

    return connections_to_atom
                    
def count_nn(nn,xyz_file_read):
    countC = 0
    countH = 0
    countO = 0
    countN = 0
    countF = 0
    for j in range(len(nn)):
        for line_index, each_line in enumerate(xyz_file_read.split('\n')):
            if line_index == int(nn[j])+1:
                if '1' in each_line[0:1]:
                    countH = countH+1
                if '6' in each_line[0:1]:
                    countC = countC+1
                if '7' in each_line[0:1]:
                    countN = countN+1
                if '8' in each_line[0:1]:
                    countO = countO+1
                if '9' in each_line[0:1]:
                    countF = countF + 1
    total = countC + countH + countO + countN + countF
    
    return total, countC, countH, countO, countN, countF

def check_element(neighbor,xyz_file_read):

    for line_index, each_line in enumerate(xyz_file_read.split('\n')):
        if line_index == neighbor + 1:
            if '1' in each_line[0:1]:
                element_present = 'H'
            if '6' in each_line[0:1]:
                element_present = 'C'
            if '7' in each_line[0:1]:
                element_present = 'N'
            if '8' in each_line[0:1]:
                element_present = 'O'
            if '9' in each_line[0:1]:
                element_present = 'F'

    return element_present