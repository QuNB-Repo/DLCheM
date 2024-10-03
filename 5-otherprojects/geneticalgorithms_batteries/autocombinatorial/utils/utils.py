import random
import fileinput
from re import L

def check_valency(position_to_add,mol_string,double):
    
    valency = 0 

    #default assumed the position is not available
    available = False

    #first check if this position is a nitrogen! 
    nitrogen_present = False
    count_lines = 0
    for line in mol_string.split('\n'):
        if count_lines == 3 + position_to_add:
            if ' N ' in line:
                nitrogen_present = True
        count_lines = count_lines + 1

    #if nitrogen, the valency counting starts at 1 to only allow up to 3 bonds (up to max of 4)
    if nitrogen_present == True:
        valency =  valency + 1
    
    if position_to_add == 6 and double == True:
        valency = valency + 2

    if position_to_add == 6 and double == False:
        valency = valency + 1

    #add all the bond numbers for this position
    count_lines = 0
    for line in mol_string.split('\n'):
        if count_lines > 3:
            if ' ' + str(position_to_add) + ' ' in line[0:7]:
                valency = valency + int(line[8])
        count_lines = count_lines + 1

    #if position has valency less than 4 then it becomes available, else stays unavailable
    if valency < 4:
        available = True
    if valency >= 4:
        available = False
    
    return available

def build_fg(mol_string, fg_data_file, position_for_fg,benchmark_atom_number,benchmark_bond_number):

    number_added_atoms = 0
    number_added_bonds = 0
    modified_string = ''
    print('before',mol_string)
    for line in mol_string.split('\n'):
        if ' 1  0 ' in line or ' 2  0 ' in line or ' 3  0 ' in line:
            for line2 in fileinput.FileInput(fg_data_file,inplace=0):
                if 'a' in line2:
                    if benchmark_atom_number >= 9:
                        line2 = line2.replace('a',' ' +str(position_for_fg))
                    elif benchmark_atom_number < 9:
                        line2 = line2.replace('a', ' ' + str(position_for_fg))
                if 'm' in line2:
                    if benchmark_atom_number >= 9:
                        line2 = line2.replace('m', str(benchmark_atom_number+1))
                    elif benchmark_atom_number < 9:
                        line2 = line2.replace('m', ' ' + str(benchmark_atom_number+1))
                if 'n' in line2:
                    if benchmark_atom_number >= 8:
                        line2 = line2.replace('n', str(benchmark_atom_number+2))
                    elif benchmark_atom_number < 8:
                        line2 = line2.replace('n', ' ' + str(benchmark_atom_number+2))
                if 'o' in line2:
                    if benchmark_atom_number >= 7:
                        line2 = line2.replace('o', str(benchmark_atom_number+3))
                    elif benchmark_atom_number < 7:
                        line2 = line2.replace('o', ' ' + str(benchmark_atom_number+3))
                if 'p' in line2:
                    if benchmark_atom_number >= 6:
                        line2 = line2.replace('p', str(benchmark_atom_number+4))
                    elif benchmark_atom_number < 6:
                        line2 = line2.replace('p', ' ' + str(benchmark_atom_number+4))
                if 'q' in line2:
                    if benchmark_atom_number >= 5:
                        line2 = line2.replace('q', str(benchmark_atom_number+5))
                    elif benchmark_atom_number < 5:
                        line2 = line2.replace('q', ' ' + str(benchmark_atom_number+5))
                if 'r' in line2:
                    if benchmark_atom_number >= 4:
                        line2 = line2.replace('r', str(benchmark_atom_number+6))
                    elif benchmark_atom_number < 4:
                        line2 = line2.replace('r', ' ' + str(benchmark_atom_number+6))
                if 's' in line2:
                    if benchmark_atom_number >= 3:
                        line2 = line2.replace('s', str(benchmark_atom_number+7))
                    elif benchmark_atom_number < 3:
                        line2 = line2.replace('s', ' ' + str(benchmark_atom_number+7))
                if ' 1  0 ' in line2 or ' 2  0 ' in line2 or ' 3  0 ' in line2:
                    number_added_bonds = number_added_bonds + 1
                if ' C  ' in line2 or ' N  ' in line2 or ' Cl  ' in line2 or ' O  ' in line2 or ' F  ' in line2 or ' S  ' in line2:
                    number_added_atoms = number_added_atoms + 1
                modified_string = modified_string + line2 
            break

        modified_string = modified_string + line + '\n' 
        
    for line in mol_string.split('\n'):
        if ' 1  0 ' in line or ' 2  0 ' in line or ' 3  0 ' in line:
            modified_string = modified_string + line + '\n'
        if 'END' in line:
            modified_string = modified_string + line + '\n'

    save_fg_index = benchmark_atom_number + 1
    benchmark_atom_number = benchmark_atom_number + number_added_atoms
    benchmark_bond_number = benchmark_bond_number + number_added_bonds

    mol_string = modified_string
    print('AFTER',mol_string)
    return mol_string, benchmark_atom_number, benchmark_bond_number, save_fg_index



def count_fg_atoms(filename):

    fg_numb_atoms = 0
    for line in fileinput.FileInput(filename,inplace=0):
        if ' C  ' in line or ' N  ' in line or ' O  ' in line or ' Cl  ' in line or ' S ' in line: 
            fg_numb_atoms = fg_numb_atoms + 1
    
    return fg_numb_atoms

def adjust(string,index,shift,number_atoms,number_bonds,shift1):



    number_atoms = number_atoms - shift
    number_bonds = number_bonds - shift1

    new_string = ''
    count_lines = 0 
    for line in string.split('\n'):
        full_line = line
        if (number_atoms+number_bonds+3) >= count_lines > (number_atoms+3) and 'END' not in line:
            line1 = line.replace(line,line[0:3])
            line2 = line.replace(line,line[3:12])
            if int(line1[1:3]) > index:
                if int(line1[1:3])-shift1 < 10:
                    line1 = line1.replace(line1[1:3],' '+str(int(line1[1:3])-shift1))
                if int(line1[1:3])-shift1 >= 10:
                    line1 = line1.replace(line1[1:3],str(int(line1[1:3])-shift1))
            if int(line2[1:3]) > index:
                if int(line2[1:3])-shift1 < 10:
                    line2 = line2.replace(line2[1:3],' '+str(int(line2[1:3])-shift1))
                if int(line2[1:3])-shift1 >= 10:
                    line2 = line2.replace(line2[1:3],str(int(line2[1:3])-shift1))
            full_line = line1 + line2
        new_string = new_string + full_line + '\n'
        count_lines = count_lines + 1
 
    final_string = ''
    count_lines = 0
    for line in new_string.split('\n'):
        if '0999' in line:
            line = line.replace(' '+ line[1:3]+' '+line[4:6],' ' + str(number_atoms)+' ' + str(number_bonds))
        final_string = final_string+ line + '\n'
        count_lines = count_lines + 1
    
    return final_string
        

                
