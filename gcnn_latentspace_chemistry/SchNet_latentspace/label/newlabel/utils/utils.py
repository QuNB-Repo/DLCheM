import os

import random

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

 #   return mol_molecule

def labelmolfile(elements,expandinter,fg_dictionary,fg_label_list,fg_labellda_dictionary,fg_labellda_list,count_fg_index,fg_marker_dictionary,fg_marker_list):

    #define random decimal color generator
    def random_decimal_color():
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)
        decimal_color = (r << 16) + (g << 8) + b
        return decimal_color

    #define the label function, returns a fg key for each fg found in the atomistic dataset
    def label_fg(atom_index,number_atoms,atom_id):
        #go through the bonding section of mol file, find neighbor atoms,
        #go through a copy of the mol file, find the neighbor atom's identity
        #write a string of ALL the neighboring atom's identities and order string alphabetically
        #put string in dictionary as a key with the value as a decimal color for gnuplot


        if expandinter == True:
            for line_index, line in enumerate(mol_string.split('\n')):
                if line_index == atom_index + 3:
                    atom_id = line[31:32]
            fg_string = atom_id + '-'
        else:
            fg_string = atom_id+'-'

        
        unordered_subfg_string = ''
        for line_index, line in enumerate(mol_string.split('\n')):
            
            #in bonding section of mol file
            if line_index > 3 + number_atoms:

                #if atom_index in line
                if ' '+str(atom_index)+' ' in line[0:7]:
                    if ' '+str(atom_index) +' ' in line[0:4]:
                        neighbor_index = int(line[3:7])
                    if ' '+str(atom_index)+' ' in line[3:7]:
                        neighbor_index = int(line[0:4])
            
                    
                    for line_index_of_copy, line_of_copy in enumerate(mol_string.split('\n')):
                        
                        if line_index_of_copy == 3 + neighbor_index:
                            

                            element = line_of_copy[31:32]
                    
                            unordered_subfg_string = unordered_subfg_string + element
                            ordered_subfg_string = ''.join(sorted(unordered_subfg_string))
        
        fg_string = fg_string + ordered_subfg_string

#        fg_string = fg_string[0:2].join(sorted(fg_string[2:len(fg_string)-1]))

        #run deeper label again for certain groups
    #    if fg_string == 

        return fg_string


    #ALLOWED marker types
    marker_types = [1,7,8,10,13,2,5,11,20,12,15]

    elements_dictionary = {}
    for element in elements:
        if element == 1:
            key = 'H'
            value = 1
            elements_dictionary[key] = value
        if element == 6:
            key = 'C'
            value = 6
            elements_dictionary[key] = value
        if element == 7:
            key = 'N'
            value = 7
            elements_dictionary[key] = value
        if element == 8:
            key = 'O'
            value = 8
            elements_dictionary[key] = value
        if element == 9:
            key = 'F'
            value = 9
            elements_dictionary[key] = value


    #open temporary mol file
    mol_file = open('temp.mol',mode='r')

    #read into a string
    mol_string = mol_file.read()


            
        
    #parse into lines and look for the elements that need to be analyzed
    for line_index, line in enumerate(mol_string.split('\n')):
        if line_index == 3:
            number_atoms = int(line[1:3])

        if line_index > 3:
            if line[31:32] in elements_dictionary:

                #label that atom's functional group, if it is a new connection add it to a fg_dictionary 
                #name starting from central atom
                atom_index = line_index - 3
                atom_id = line[31:32]

                #go through the mol_string again to find the connectivity of the atom and the elements
                #if dictionary_element is new add to a dictionary set



                #if you need the interacting atom's functional froups then run label_fg on each of the interacting atom_indices
                if expandinter == True:
#                    for line_index_of_copy, line_of_copy in enumerate(mol_string.split('\n')):
             
                            #if atom_index in line
                            #ACTUALLY ALL OF THEM ARE NEIGHBORS!! except for itself

                            #use this if you want only direct neighbors
#                        if line_index_of_copy > 3 + number_atoms:                           
#                            if ' '+str(atom_index)+' ' in line_of_copy[0:7]:
#                                print('line_of_atom_index',line_of_copy)
#                                if ' '+str(atom_index)+' ' in line_of_copy[0:4]:
#                                    neighbor_index = int(line_of_copy[3:7])
#                                if ' '+str(atom_index)+' ' in line_of_copy[3:7]:
#                                    neighbor_index = int(line_of_copy[0:4])

                    for each_neighbor in range(1,number_atoms+1):
                        if each_neighbor != atom_index:
                            neighbor_index = each_neighbor
                            fg_string = label_fg(neighbor_index,number_atoms,atom_id)
                            
                            random_color = random_decimal_color()
                            if fg_string not in fg_dictionary.keys():

                                fg_dictionary[fg_string] = random_color

                                count_fg_index = count_fg_index + 1
                                fg_labellda_dictionary[fg_string] = count_fg_index

                                fg_marker_dictionary[fg_string] = random.choice(marker_types)

                            fg_label_list.append(fg_dictionary[fg_string])
                            fg_labellda_list.append(fg_labellda_dictionary[fg_string])
                            fg_marker_list.append(fg_marker_dictionary[fg_string])
                #if not expandinter then just label the atom_index as standard
                else:
                    fg_string = label_fg(atom_index,number_atoms,atom_id)

                    random_color = random_decimal_color()
                    if fg_string not in fg_dictionary.keys():

                        fg_dictionary[fg_string] = random_color

                        count_fg_index = count_fg_index + 1
                        fg_labellda_dictionary[fg_string] = count_fg_index
                        
                        fg_marker_dictionary[fg_string] = random.choice(marker_types)


                    fg_label_list.append(fg_dictionary[fg_string])
                    fg_labellda_list.append(fg_dictionary[fg_string])

    return fg_dictionary, fg_label_list, fg_labellda_dictionary, fg_labellda_list, fg_marker_dictionary, fg_marker_list, count_fg_index