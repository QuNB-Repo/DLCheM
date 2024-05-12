import os

import fileinput

import numpy as np


#Writes xyz file from a loaded database (db), needs props of db as input, and 
#required for labelling code
def write_xyz_from_db(props,idx,file_path):
    
    name_xyz = file_path+str(idx)+'.xyz'
    # open an empty file
    xyz_file = open(name_xyz,mode='w',encoding='utf-8')

    # copy props['_atomic_numbers'], props['_positions] tensor and change tensor to numpy array 
    atomic_numbers = props['_atomic_numbers']
    atomic_numbers = atomic_numbers
    
    number_atoms = len(atomic_numbers)

    positions = props['_positions']
    
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
    
    output = os.system('obabel ' + name_xyz + ' -O ' + name_mol + '>/dev/null 2>&1')

    return name_mol

def connections_list(name_mol,number_atom):

    #initialize connections list which will hold each atom's connections
    connect_list = []

    #initialize atomwise element dictionary
    element_list = []
    #initialize atom_count
    count_atom = 1

    #count lines in file
    count_line = 0
    for line in fileinput.FileInput(name_mol,inplace=0):
        count_line = count_line + 1

        #at the coordinates section
        if  5 + number_atom > count_line >= 5:
            
            #take each element
            atomwise_element = line[31:32]

            #add in atomwise_elements list within list
            connect_list.append([count_atom])
            count_atom = count_atom+1

            #add to atomwise_elements list
            element_list.append([atomwise_element])


        if count_line >= 5 + number_atom:

            #delete parts of the line that will conflict with later reading
            line = line.replace(line,line[:7])

            for j in range(number_atom):
                if ' ' + str(j+1) + ' ' in line and 'RAD' not in line:
                    if str(line[1:4]) == ' '+str(j+1)+' ' or str(line[0:4]) ==' '+str(j+1)+' ':
                        connect_list[j].append(int(line[4:7])) 
                        index = int(line[4:7]) - 1
                        element_list[j].append(element_list[index][0])
                    else:
                        connect_list[j].append(int(line[1:4]))
                        index = int(line[1:4]) - 1
                        element_list[j].append(element_list[index][0])
    
    return connect_list, element_list

def order_branches(connect_list,element_list,n_branches):
    
    if n_branches > 1:
        
        orderedelement_list = []
        orderedconnect_list = []
        
        #get fluorines first
        for i in range(len(element_list)):
            if element_list[i][0] == 'F':
                orderedelement_list.append(['F',element_list[i][1:len(element_list)]])
                orderedconnect_list.append([i+1,connect_list[i][1:len(connect_list)]])


        #get oxygens second
        for i in range(len(element_list)):
            if element_list[i][0] == 'O':
                orderedelement_list.append(['O',element_list[i][1:len(element_list)]])
                orderedconnect_list.append([i+1,connect_list[i][1:len(connect_list)]])

        #get nitrogens third
        for i in range(len(element_list)):
            if element_list[i][0] == 'N':
                orderedelement_list.append(['N',element_list[i][1:len(element_list)]])
                orderedconnect_list.append([i+1,connect_list[i][1:len(connect_list)]])

        #get carbons
        for i in range(len(element_list)):
            if element_list[i][0] == 'C':
                orderedelement_list.append(['C',element_list[i][1:len(element_list)]])
                orderedconnect_list.append([i+1,connect_list[i][1:len(connect_list)]])

        #get hydrogens
        for i in range(len(element_list)):
            if element_list[i][0] == 'H':
                orderedelement_list.append(['H',element_list[i][1:len(element_list)]])
                orderedconnect_list.append([i+1,connect_list[i][1:len(connect_list)]])

        #now order all the atoms in each branch
        for j in range(n_branches):
        
            neighbor = orderedelement_list[j][1]

            ordered_neighborconnect = []
            ordered_neighborelement = []


            #get fluorines first
            for i in range(len(neighbor)):
                if neighbor[i] == 'F':
                    ordered_neighborelement.append('F')
                    ordered_neighborconnect.append(orderedconnect_list[j][1][i])
            #get oxygens second
            for i in range(len(neighbor)):
                if neighbor[i] == 'O':
                    ordered_neighborelement.append('O')
                    ordered_neighborconnect.append(orderedconnect_list[j][1][i])

            for i in range(len(neighbor)):
                if neighbor[i] == 'N':
                    ordered_neighborelement.append('N')
                    ordered_neighborconnect.append(orderedconnect_list[j][1][i])

            for i in range(len(neighbor)):
                if neighbor[i] == 'C':
                    ordered_neighborelement.append('C')
                    ordered_neighborconnect.append(orderedconnect_list[j][1][i])

            for i in range(len(neighbor)):
                if neighbor[i] == 'H':
                    ordered_neighborelement.append('H')
                    ordered_neighborconnect.append(orderedconnect_list[j][1][i])

            #replace ordered for the one in the original list
            orderedelement_list[j][1] = ordered_neighborelement
            orderedconnect_list[j][1] = ordered_neighborconnect
    
    
    if n_branches == 1:

        orderedconnect_list = connect_list
        orderedelement_list = element_list

        neighbor = orderedelement_list

        ordered_neighborconnect = []
        ordered_neighborelement = []

        #get fluorines first
        for i in range(len(neighbor)):
            if neighbor[i] == 'F':
                ordered_neighborelement.append('F')
                ordered_neighborconnect.append(orderedconnect_list[i])
        #get oxygens second
        for i in range(len(neighbor)):
            if neighbor[i] == 'O':
                ordered_neighborelement.append('O')
                ordered_neighborconnect.append(orderedconnect_list[i])

        for i in range(len(neighbor)):
            if neighbor[i] == 'N':
                ordered_neighborelement.append('N')
                ordered_neighborconnect.append(orderedconnect_list[i])

        for i in range(len(neighbor)):
            if neighbor[i] == 'C':
                ordered_neighborelement.append('C')
                ordered_neighborconnect.append(orderedconnect_list[i])

        for i in range(len(neighbor)):
            if neighbor[i] == 'H':
                ordered_neighborelement.append('H')
                ordered_neighborconnect.append(orderedconnect_list[i])
        
        orderedelement_list = ordered_neighborelement
        orderedconnect_list = ordered_neighborconnect

    return orderedconnect_list, orderedelement_list
        #get oxygens first
#        for i in range(len(connecting_atoms)):
#            if connecting_atoms[i] == 'O':
#                ordered_connecting_atoms.append('O')

def perturbation(totalconnect_list,totalelement_list,connect_list,element_list,n_branch,save_root,n_pert):

    fullconnect_list = []
    fullelement_list = [] 
    fullconnect_list.append([connect_list])
    fullelement_list.append([element_list])

    if n_pert == 0:
        
        #define root atom as the fullconnect_list
        fullconnect_list = [[[connect_list[0]]]]
        fullelement_list = [[[element_list[0]]]]


    if n_pert == 1:
        #define root atom as the fullconnect_list
        fullconnect_list = [[connect_list]]
        fullelement_list = [[element_list]]


    else:
        for pert in range(n_pert):
            newconnect_list = []
            newelement_list = []

            if n_branch == 1:

                #initiate neighbor list
                neighborconnect_list = []
                neighborelement_list = []

                if pert > 0:
                    n_neighbor = len(connect_list[0][0][1])
                else:
                    n_neighbor = len(connect_list[1])


                for neighbor in range(n_neighbor):
                    
                    if pert > 0:
                        nb = connect_list[0][0][1][neighbor]
                    else:
                        nb = connect_list[1][neighbor]


                    if nb not in save_root:
                        save_root.append(nb)

                        #find each neighbor's connect_list and element_list and store them in a set of lists
                        for connection in range(len(totalconnect_list)):
                            if totalconnect_list[connection][0] == nb:
                                neighborconnect_list.append(totalconnect_list[connection])
                                neighborelement_list.append(totalelement_list[connection])

                newconnect_list.append(neighborconnect_list)
                newelement_list.append(neighborelement_list)

                #append neighbor's connections to new_connect list
                fullconnect_list.append(neighborconnect_list)
                fullelement_list.append(neighborelement_list)

                #newconnect_list now becomes the "old" conncet_list and we try again for the next pert
                connect_list = newconnect_list 
                element_list = newelement_list
                n_branch = len(connect_list[0])
                #first index is number of perts done, second is number of branches in pert, third is root vs neighbor
            elif n_branch > 1:
                neighborconnect_list = []
                neighborelement_list = []
                #for each branch, do the above, but with slight twist, not allowed to go back! save roots that you have already done
                for branch in range(n_branch):

                    
                    #define branch
                    brch = connect_list[0][branch]
                    
                    #define neighbor list
                    nbr_list = connect_list[0][branch][1]
                    n_neighbor = len(nbr_list)

                    for neighbor in range(n_neighbor):
                        
                        nbr = connect_list[0][branch][1][neighbor]

                        #if neighbor is not in saved roots

                        if nbr not in save_root:
                            save_root.append(nbr)
                            #find each neighbor's connect_list and element_list and store them in a set of lists
                            for connection in range(len(totalconnect_list)):
                                if totalconnect_list[connection][0] == nbr:
                                    neighborconnect_list.append(totalconnect_list[connection])
                                    neighborelement_list.append(totalelement_list[connection])
                        
                newconnect_list.append(neighborconnect_list)
                newelement_list.append(neighborelement_list)

                #append neighbor's connections to new_connect list
                fullconnect_list.append(neighborconnect_list)
                fullelement_list.append(neighborelement_list)

                #newconnect_list now becomes the "old" conncet_list and we try again for the next pert
                connect_list = newconnect_list 
                element_list = newelement_list
                n_branch = len(connect_list[0])

    return fullconnect_list, fullelement_list


def fg_key(new_connect,new_element,save_root,pert):

    if pert == 0:
        fg_key = new_element[0][0][0]
    #just order the branches according to most oxygen,nitrogen,hydrogen
    else:
        for pert in range(len(new_element)):

            #define #number of branches
            n_branch = len(new_element[pert])
            brch = new_element[pert]

            #if there are more than one branch, 
            #give a rating to each branch and reorder the branches according to rating
            rating = []
            newbranch_order = []
            if n_branch == 1:
                newbranch_order = new_element[pert]
            if n_branch > 1:
                
                for branch in range(n_branch):
                #for each neighbor list sum priority numbers
                    n_neighbor = len(new_element[pert][branch][1])

                    score = 0
                    for neighbor in range(n_neighbor):

                        if new_element[pert][branch][1][neighbor] == 'F':
                            score = score + 1000000 + 1/(branch+1)
                        if new_element[pert][branch][1][neighbor] == 'O':
                            score = score + 100000 + 1/(branch+1)
                        if new_element[pert][branch][1][neighbor] == 'N':
                            score = score + 10000 + 1/(branch+1)
                        if new_element[pert][branch][1][neighbor] == 'C':
                            score = score + 1000 + 1/(branch+1)
                        if new_element[pert][branch][1][neighbor] == 'H':
                            score = score + 100 + 1/(branch+1)
                    rating.append(score)
                
                    #go through max rating index as index for ordering each branch
                index_done = []
                rating2 = rating.copy()
                for branch in range(n_branch):
                    max_value = max(rating2)

                    max_index = rating.index(max_value)


                    count_max = rating.count(max)
                    newbranch_order.append(brch[max_index])
                    
                    max_index2 = rating2.index(max_value)
                    rating2.pop(max_index2)
            #replace pert
            new_element[pert] = newbranch_order

            #construct fg key directly from organized list
        fg_key = ''
        for pert in range(len(new_element)):

            n_branch = len(new_element[pert])
            brch = new_element[pert]

            for branch in range(n_branch):

                fg_key = fg_key + brch[branch][0] + '-'

                n_neighbor = len(new_element[pert][branch][1])
                
                fg_key = fg_key + '['

                for neighbor in range(n_neighbor):

                    fg_key = fg_key + new_element[pert][branch][1][neighbor] + '-'
                fg_key = fg_key + ']-'
    return fg_key

