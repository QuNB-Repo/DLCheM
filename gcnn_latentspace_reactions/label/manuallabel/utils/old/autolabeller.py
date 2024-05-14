'''
This code automatically constructs the functional group label dictionary
of up to arbitrary number of neighbors included in the functional group environment
(i.e looks recursively deeper into the neighbors through the molecule)
'''

import fileinput
import numpy as np

from schnetpack.datasets import QM9

from numpy import savetxt

from utils import utils

import os


#initialize label parameters
class label():
    def __init__(self,qm9data,label_filepath,start,end,element):
        self.label_filepath = label_filepath
        self.start = start
        self.end = end
        self.qm9data = qm9data
        element = ' ' + element + '  '
        self.element = element
        
        #depending on target element there is a slightly different labelling code
        if element == ' H  ':
            label_h(self)
        if element == ' C  ':
            label_c(self)
        if element == ' N  ':
            label_n(self)
        if element == ' O  ':
            label_o(self)

def label_h(self):
    start = self.start
    end = self.end
    label_filepath =self.label_filepath
    qm9data = self.qm9data
    element = self.element

    #initialize functional group data list 
    datah = list()
    #initialize output file
    output_file = open(label_filepath,mode='a',encoding='utf-8-sig')
    #initialize starting label_id
    label = 0
    #initialize functional group dictionary
    env_dict = dict()
    #set perturbation value (how many neighbors to look into successively)
    pert = 3

    #for each molecule
    for idx in range(start,end):

        #get molecule properties from qm9 dataset
        at, props = qm9data.get_properties(idx)
        #get number of atoms
        number_atoms = len(props['_atomic_numbers'])
        #get atomic numbers
        atomic_number = props['_atomic_numbers']

        #write xyz from db
        label_direc = './'
        name_xyz = utils.write_xyz_from_db(props,idx,label_direc)
        #convert xyz to mol
        name_mol = utils.xyz_to_mol(idx,label_direc)

        #check that H exists in the molecule
        if 1 in props['_atomic_numbers']:

            #get hydrogen indices
            all_hydrogen_indices = utils.store_positions(idx,name_mol,props,element)
            
            #get their immediate connections (connections matrix)
            all_element_connections = utils.connection_matrix(all_hydrogen_indices,name_mol,number_atoms)
            print('all_H', all_hydrogen_indices)
            print('all_connected_to_H',all_element_connections)
            
            #start with each element in the molecule as "root" for the functional_group key
            #for each element's connnections
            for j in range(len(all_element_connections)):

                #define this atom's connections
                this_atom_connections = all_element_connections[j]

                #get connecting elements counts
                total, countC, countH, countO, countN, countF = utils.count_nn(this_atom_connections,name_xyz)
                
                #count all non-H's, determine number of branches 
                this_atom_n_branches = total

                #effective_branches
                effective_branches = countC + countN + countO
                
                #construct fg key from counts
                funcgroup_key = utils.create_fglabel(total,countH,countC,countN,countO,countF,effective_branches)

                #run each pert for each branch, adding a subscript for branch number
                save_neighbor_indices = [all_hydrogen_indices[j]]

                #begin perturbation count
                count_pert = 1

                #figure out which atoms exist and on which index
                elements_present = []

                for b in range(this_atom_n_branches): 
 
                    #define the neighbor
                    neighbor = this_atom_connections[b]
                    
                    #check what the element the neighbor is
                    element = utils.check_element(neighbor,name_xyz)

                    #append element to neighbor element list
                    elements_present.append(element)

                #find which of the branches is oxygen
                for b in range(this_atom_n_branches):

                    #define the neighbor
                    neighbor = this_atom_connections[b]

                    #check that the neighbor/branch has not been already evaluated
                    if neighbor not in save_neighbor_indices:
                    
                        #first do the oxygen perturbations all in a bracket (diff syntax to funcgroup_key)
                        if elements_present[b] == 'O':
                            
                            n_pert == 1 
                            while n_pert < pert: 

                                #run perturbation code (runs exactly same code, different name, but now on this atom)
                                funcgroup_key = utils.pert_code()

                                n_pert = n_pert + 1

                        

                #check if not in dictionary
                if funcgroup_key not in env_dict:
                    label = label + 1 
                    env_dict[funcgroup_key] = label

                if funcgroup_key in env_dict:
                    label = env_dict[funcgroup_key]

                row = '' +str(label)+' '+str(idx)+' '+funcgroup_key+' \n'
                output_file.write(row)

    output_file.close()