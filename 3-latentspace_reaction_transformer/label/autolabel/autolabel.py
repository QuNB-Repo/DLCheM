from utils import utils_connect

from schnetpack.datasets import QM9
import numpy as np
from numpy import genfromtxt, savetxt
from schnetpack import AtomsData
'''
Last Updated: 2024-02-20

LABELING ATOMIC ENVIRONMENTS AUTOMATICALLY AT ANY DEPTH FROM AN ATOMIC CENTER

This code handles labelling atoms by their atomic environment in a molecular database, db file
It can do this up to any DEPTH of bonds specified (solves a combinatorial problem in real time!)

the output of the method is a csv file of columns:
            label | mol_index | atom_index

            output_filepath     - filepath which will contain the csv file described above
            start               - beginning index to label the dataset
            end                 - end index of labelling the dataset
            qm9                 - boolean if qm9 db is being used, load it with QM9 method from SchNet,
                                  otherwise use AtomsData
            dataset_filepath    - db file that is to be labelled
            n_depth             - depth of labelling, how many bonds in to differentiate between atomic environments
                                  depth = 1 means one bond away (H-O-, from persp. of H), 2 means two bonds (H-O-C-, from persp. of H), ... etc
            elements            - target elements to label 
            scratch_direc       - scratch directory
'''
def label(output_filepath,start,end,qm9,dataset_filepath,n_depth,elements,scratch_direc):
    '''
    env_dict            - this environment dictionary keeps track of the label of environments that have been found
                          so that if the environment is found again, it is labelled according to the dictionary,
                          otherwise a new fg_key is added to the dictionary with the new fg_label
    unique_label_count  - a counter that keeps track of how many atomic environments were found, 
                          also used to integerize the fg_key for each atomic environment founr
    output_file         - opens the output_filepath which will hold the output csv descrived above
                          in "append" mode so that it can be continued if interrupted
    each_molecule       - the currect integer of the dataset (runs through integers of the dataset)
    data                - loading the dataset dictionary using QM9 method or AtomsData method
    at, prop            - loads molecular formula and properties dictionary of the dataset
    number_atom         - number of atoms in each_molecule
    atomic_number       - the list of atomic numbers in the molecule
    utils_connect.
    write_xyz_from_db   - this executes a utility a scratch xyz file from db file to be converted to mol file
    utils_connect.
    xyz_to_mol          - this executes a utility that converts xyz file to mol file
    utils_connect.      
    connections_list    - this executes a utility that outputs the connect list and element, see below
    connect_matrix      - (target_atom,neighbor_idx) a matrix where the rows are the target_atoms successively as they appear in QM9 (for labelling)
                          and the columns represent the neighbor indices of those targets
                          NOTE that for each depth this gets replaced/updated with new targets/neighbors 
                          (the neighbors of previous depth become the targets of next depth...)
    element_matrix      - (target_atom,neighbor_element) a matrix where the rows are the target_atoms succesively as they appear in QM9 (for labelling)
                          and the columns represent the neighbor elements of those targets
                          NOTE that for each depth this gets replaced/updated with new targets/neighbors 
                          (the neighbors of previous depth become the targets of next depth...)
    utils_connect.      
    connections_list    - this executes the tool that finds the matrices above using a mol file, it does this to initiate the first depth
                          and is repeated at every depth to update the matrices using new targets (new neighbors of prev depth become the targets of next depth)
    n_targets           - this is a counter that also gets updated keeping track of the number of targets at each depth
                          THIS IS REALLY ONLY IMPORTANT WHEN YOU ARE NOT TARGETTING FOR ALL ELEMENTS IN DATASET. OTHERWISE N_BRANCHES = N_ATOMS.
                          if you are targetting only H then it is important to keep track of how many H's there are
    adj_idx_matrix     - this is the main tool used by the algorithm. EACH row of the matrix
                         will serve as a root target atom for labelling. Each row will go into the depth labelling to get labelled
    adj_ele_matrix     - this is the same as above but instead of indices it is the elements of the atoms
  
                          NOTE BOTH MATRICES BUT are ordered in a particular way so the label always comes out consistent and based on priority! 
                          that is if your connections are H,C,N,O, then you are the same as C,H,N,O (Switching)... BUT THIS only works if we order first, O,N,C,H (we chose atomic number priority)
                          WE LABEL biggest atomic numbers first
                          IT IS ORDERED according to the function below
    utils_connect.
    order_branches      - this executes the tool that ORDERS tha adhaceny matrix according to priority      
                          the order is according to priority of atomic_number (similar to how it is done in chemistry nomenclature),
                          F gets higher priority than O, H gets least priority at each depth
                          this priority system holds BOTH for the root atoms and the neighbors of the atom
                          the final "tensors" look like this:

                            [ [1, [3,2]],  [3,[1,2]], ...]
                            [['O', ['C','H']], ['C',['O','H']] ... ETC  

                          so that each "first" element of this tensor looks like ['O',['C','H']]  

                          such a tensor will be used to make the final fg_keys for each target root/element/atom/row of adjacency along its depth 
                          that will differentiate each atom's environment from the rest using environments dictionary   
    each_target         - an integer that runs through all the targets at a certain depth
    
    target_connect_list - runs through the connection tensor described above and takes each branch/target, ex. [1,[3,2]]
    target_element_list - runs through the element branches described above and takes each branch/target, ex. ['O',['C','H']]
    utils_connect.
    depth_label         - the main code that finds the depth of branches/neighbors using the adjacency
                          returns the result as target_depth_elements
    target_depth_
    elements            - the result of exploring depths of a target atom as described above
                          returns 
                          
                          [ (depths) [ (branches) [ roots, [root_connections]]]]
    
    utils_connect.
    order_branches_
    in_each_depth       - goes through each depth and orders the branches based on a atomic priority rating 
                         (a sum of exponentially increasing atomic number priority all neighbors of that branch)
    utils_connect.
    fg_key              - writes the fg_key based on ordered branches of the target_depth_element
                          always follows this format:

                         'rootelement-[branching_nbrs]-rootelement-[branching_nbrs]'

                         BUT ROOTELEMENT IS ORGANIZED ACCORDING TO DEPTH
                         AND BRANCHING_NBRS ARE ORGANIZED ACCORDING TO ATOMIC NUMBER PRORITY

                         inconsequentially, targets are also organized according to atomic number priority
    row                - string row that will go in output file labelling the target atom
    label              - the label integer of the fg from the unique_label_count
    inset              - a boolean "inside" or "outside" set to help identify whether a target 
                         was the first time encountered "out" 
                         or was already "in" the set
    '''
    #initiate chemical environments dictionary
    env_dict = dict()
    #initiate chemical environment count
    unique_label_count = 0 

    #open output file in "append" mode
    output_file = open(output_filepath,mode='a',encoding='utf-8-sig')

    #write a title to the csv file we are building
    title_row = 'label,mol_idx,atom_idx,key,in/out \n'

    output_file.write(title_row)
    
    #run through the specified start-end integers of the dataset
    for each_molecule in range(start,end):

        if each_molecule % 1000 == 0:
            print(each_molecule)

        #Load qm9data
        if qm9 == True:
            data = QM9(dataset_filepath,download=False,remove_uncharacterized=True)
        else:
            data = AtomsData(dataset_filepath,available_properties=['energy'])

        #get molecule properties from qm9 dataset
        at, prop = data.get_properties(each_molecule)
        #get number of atoms
        number_atoms = len(prop['_atomic_numbers'])

        #write xyz from db
        utils_connect.write_xyz_from_db(prop,scratch_direc)
        #convert xyz to mol
        utils_connect.xyz_to_mol(scratch_direc)

        #construct total_branches_list
        name_mol = scratch_direc + 'temp.mol'
        adj_ind_matrix, adj_ele_matrix = utils_connect.build_adjacency_matrix(name_mol,number_atoms)

        #ad ind matrix
        n_targets=len(adj_ind_matrix)
        adj_ind_matrix, adj_ele_matrix = utils_connect.order_adjacency(adj_ind_matrix, adj_ele_matrix,n_targets)
        

        '''
        THE CODE IS BASED ON THE CONCEPT OF 
        "INITIALIZING A TARGET ROOT/ATOM's LABEL"              - for each atom in the molecule, we define a root branch, WHICH IS JUST A ROW IN THE ADJACENCY
        "LABEL THE DEPTHS OF EACH ROOT/ATOM/ROW of ADJACENCY"  - each target atom/element will be labelled --> so for each:
                                                                 run through the neighbors of the root atom/element, and find their neighbors' neighbors, 
                                                                 and keep going (depths), using the adjacency to explore
                                                                 and use this to LABEL EACH ATOM". All is ordered according to priority
        
        NOTE it cannot be initialized where all target elements are done at once. 
        It needs to root each target and treat separately
        '''

        for each_target in range(n_targets):
            #define root branch
            target_connect_list = adj_ind_matrix[each_target]
            target_element_list = adj_ele_matrix[each_target]

            #run perturbation labelling on target through the depth specified to get
            #to get the target_depth_connections, target_depth_elements, all will be used to make the fg_key
            target_depth_elements = utils_connect.depth_label(adj_ind_matrix,adj_ele_matrix,target_connect_list,target_element_list,n_depth)

            #order branches of each depth in the target_depth_elements
            target_depth_elements = utils_connect.order_branches_in_each_depth(target_depth_elements,n_depth)


            #make the fg_key using the target_depth_element fully ordered
            key = utils_connect.fg_key(target_depth_elements,n_depth)


            if key in env_dict:
                if str(key[0]) in elements:
                    label = env_dict[key]
                    inset = 'in'

            if key not in env_dict:
                if str(key[0]) in elements:
                    unique_label_count = unique_label_count + 1 
                    label = unique_label_count
                    env_dict.update({key: unique_label_count})
                    inset = 'out'
            
            #change this to allow any element(s) possibility
            if str(key[0]) in elements:
                row = '' +str(label)+','+str(each_molecule)+','+str(each_target)+','+key+','+ str(inset) +  '\n'
                output_file.write(row)

    output_file.close()

    return unique_label_count
