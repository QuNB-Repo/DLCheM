import os
import fileinput
import numpy as np
'''
Last Updated: 2024-02-24

Utility functions that make the autolabel code happen! 
(labeling atomic environment automatically at any depth from atomic center)

    write_xyz_from_db              - exctracts molecules from db file as a temporary xyz file in a scratch dir
    xyz_to_mol                     - converts the temporary xyz file to a mol file 
                                     so that we can build adjacency matrix
    build_adjacency_matrix         - builds adjacency matrix of each atom in the molecule:
                                     [[1,...(bonded neighbors)],[2,...],[3,...],...,[N,...]]
    order_adjacency                - orders adjacency and reformats the data structure to better suit the problem
                                     1) priority ordering        - all roots are ordered 
                                     2) format data structure    - the data becomes [(branch)[root,[root_connections]],...,]
    depth_label                    - the main chunk of the code that derives the essential target_depth_elements
                                     these are the elements of the branches from the targer ordered according to depth/priority
    order_branches_in_each_depth   - A second round of ordering that ensures that all branches at each depth 
                                     are ordered according to score rating based on summing the priority of neighbors in each branch
                                     this ensures there is a consistent way of ordering branches so that two equivalent molecules are counted as different because of ordering
    fg_key                         - the ultimate label for the atomic environment 
                                     uses target_depth_elements with ordered branches in each depth
                                     
                                     fg_key always follows this format:

                                        'rootelement-[branching_nbrs]-rootelement-[branching_nbrs]'
'''


#Writes xyz file from a loaded database (db), needs props of db as input, and 
#required for labelling code
def write_xyz_from_db(props,scratch_direc):
    '''
    utility that helps to extract xyz information in a temporary file from db, 
    this is just so we can convert xyz to mol file, mol file help us with 
    bonding and connectivity (so we can find bonded neighbors), making the adjacency matrix

        name_xyz            - name of scratch xyz file (uses scrach_direc)
        xyz_file            - opening the xyz file in write mode
        atomic_numbers      - atomic numbers of the xyz from props dictionary 
                              (detach to numpy on earlier schnetpack == 0.3)
        number_atoms        - number of atoms from atomic numbers 
        positions           - positions of xyz from props dictionary 
                              (detach to numpy on earlier schnetpack == 0.3)
        each_atom           - integer running through each atom, checking atomic number
                              and then writing the atom's row in the xyz file labelled by element 

    
    '''
    name_xyz = scratch_direc+'temp.xyz'
    # open an empty file
    xyz_file = open(name_xyz,mode='w',encoding='utf-8')

    # copy props['_atomic_numbers'], props['_positions] tensor and change tensor to numpy array 
    atomic_numbers = props['_atomic_numbers']
    atomic_numbers = atomic_numbers.detach().numpy()
    
    number_atoms = len(atomic_numbers)

    positions = props['_positions'].detach().numpy()
    
    # write xyz file in xyz file format
    xyz_file.write(str(number_atoms)+'\n')
    xyz_file.write('Title'+'\n')
    for each_atom in range(number_atoms):
        if atomic_numbers[each_atom] == 1:
            xyz_file.write('H ' + str(positions[each_atom][0]) + ' ' + str(positions[each_atom][1]) + ' ' + str(positions[each_atom][2]) + '\n')
        if atomic_numbers[each_atom] == 6:
            xyz_file.write('C ' + str(positions[each_atom][0]) + ' ' + str(positions[each_atom][1]) + ' ' + str(positions[each_atom][2]) + '\n')
        if atomic_numbers[each_atom] == 7:
            xyz_file.write('N ' + str(positions[each_atom][0]) + ' ' + str(positions[each_atom][1]) + ' ' + str(positions[each_atom][2]) + '\n')
        if atomic_numbers[each_atom] == 8:
            xyz_file.write('O ' + str(positions[each_atom][0]) + ' ' + str(positions[each_atom][1]) + ' ' + str(positions[each_atom][2]) + '\n')
        if atomic_numbers[each_atom] == 9:
            xyz_file.write('F ' + str(positions[each_atom][0]) + ' ' + str(positions[each_atom][1]) + ' ' + str(positions[each_atom][2]) + '\n')
    xyz_file.close()


#runs obabel conversion tool on the created xyz file to convert to mol file
#required for labelling code
def xyz_to_mol(scratch_dir):
    '''
    utility that converts temporary xyz file to temporary mol file
    because mol file helps us constrcut the adjacency matrix

        name_xyz        - name of temporary xyz file (needs scratch_dir)
        name_mol        - name of the temporary mol_file (needs scratch_dir)
        output          - a variable to print to ensure that the os operation is working, 
                          if it prints 1 then it worked, if 0 then it did not work
    '''
    # define name of temporary xyz file according to idx, define the name of temporary mol file 
    name_xyz = scratch_dir + 'temp.xyz'
    name_mol = scratch_dir  + 'temp.mol'
    
    #use obabel to convert xyz to mol
#    output = subprocess.run('obabel ' + name_xyz + ' -O ' + name_mol)

    output = os.system('obabel ' + name_xyz + ' -O ' + name_mol )# + ' >/dev/null 2>&1')

    return name_mol

def build_adjacency_matrix(name_mol,number_atom):
    '''
    this uses the mol file to build the adjacency matrix
    the adjacency matrix is essential to exploring neighbors and the connections 
    to those neighbors in the depth-labeling code

        adj_idx_matrix        - this matrix holds the adjacency indices 
        adj_ele_matrix        - this matrix holds teh adjacency element symbols
        count_atom            - initial atom count is 1 (arbitrary, can be zero),
                                keeps track of atom count through reading the mol file
        count_line            - initial line count is 1 (arbitrary)
                                keeps track of line count through reading the mol file 

                                once count_line == 5,                     - beginning of the coordinates section
                                once 5 + number_atoms and NOT including   - end of coordinates section

                                once 5 + number_atoms                     - beginning of the bonding section
                                once 5 + number_atoms + number_bonds      - end of bonding section 

        line                  - a line of the string input of molfile
        atomwise_element      - the part of the line in the coordinates section that reveals the element
        each_atom             - integer that runs through the atom indices through to number of atoms
        nbr_idx               - the part of the line in the bonding section that reveals 
                                the neighbor index this atom_idx is connected to
    '''
    #initialize connections list which will hold each atom's connections
    adj_idx_matrix = []

    #initialize atomwise element dictionary
    adj_ele_matrix = []
    #initialize atom_count
    count_atom = 1
    #count lines in file
    count_line = 1
    for line in fileinput.FileInput(name_mol,inplace=0):

        #at the coordinates section
        if  5 + number_atom > count_line >= 5:
            
            #take each element
            atomwise_element = line[31:32]

            #add in atomwise_elements list within list
            adj_idx_matrix.append([count_atom])
            count_atom = count_atom+1

            #add to atomwise_elements list
            adj_ele_matrix.append([atomwise_element])


        if count_line >= 5 + number_atom:

            #delete parts of the line that will conflict with later reading
            line = line.replace(line,line[:7])

            for each_atom in range(number_atom):
                if ' ' + str(each_atom+1) + ' ' in line and 'RAD' not in line:
                    if str(line[1:4]) == ' '+str(each_atom+1)+' ' or str(line[0:4]) ==' '+str(each_atom+1)+' ':
                        adj_idx_matrix[each_atom].append(int(line[4:7])) 
                        nbr_idx = int(line[4:7]) - 1
                        adj_ele_matrix[each_atom].append(adj_ele_matrix[nbr_idx][0])
                    else:
                        adj_idx_matrix[each_atom].append(int(line[1:4]))
                        nbr_idx = int(line[1:4]) - 1
                        adj_ele_matrix[each_atom].append(adj_ele_matrix[nbr_idx][0])
        
        count_line = count_line + 1

    return adj_idx_matrix, adj_ele_matrix

def order_adjacency(adj_idx_matrix,adj_ele_matrix,n_targets):
    '''
    Does an initial ordering of the adjacency matrix 
    (NOTE there is another ordering of branches of label after, based on score, this is NOT the same)
    The ordering is to ensure that we order higher elements first (not that important), 
    more importantly the neighbors of each atom are also ordered according to priority 
    making a unique labelling for each set of neighbors


        n_targets              - number of targets/atoms in the adjacency 
                                 (this can be maximum dimension of adjacency or less depending on what elements are chosen to label)
        ord_adj_ele_matrix     - the matrix that will hold all the elements of adjacency ORDERED
                                 the format of the matrix has also changed:
                                    
                                    [(target_atom/branch)[(root),[root_connections]],[,[]], ... , [,[]]]

        ord_adj_idx_matrix     - the matrix that will hold all the indices of the adjacency ORDERED
                                 the format of the matrix has also changed:
                                
                                    [(target_atom/branch/target)[(root),[root_connections]],[,[]], ... , [,[]]]

                                    A format that better streamlines the construction of the fg_key

        each_atom              - an integer running through all atoms of the adjacency, 
                                 used to access the adjacency
        neighbors              - the list neighbors elements for each atom/target as defined by the ordered format described above
                                 this list will be used to check which element is the neighbor
                                 and which should go first, it will check for 'F' first to fill 
                                 the ord_nbr_ele_list and ord_nbridx_list below, and then 'O', and then... highest to lowest
        ord_nbridx_list        - the list that will hold the neighbor idxs of a target/atom ordered
        ord_nbrele_list        - the list that will hold the neighbor elements of a target/atom ordered
                                 this list replaces the list of neighbors in the ordered format [each_atom/target][1], [(each_atom/branch/target)[ atom ,[neighbors]]]
    '''
    if n_targets > 1:
        
        ord_adj_ele_matrix = []
        ord_adj_idx_matrix = []
        
        #get fluorines first
        for each_atom in range(len(adj_ele_matrix)):
            if adj_ele_matrix[each_atom][0] == 'F':
                ord_adj_ele_matrix.append(['F',adj_ele_matrix[each_atom][1:len(adj_ele_matrix)]])
                ord_adj_idx_matrix.append([each_atom+1,adj_idx_matrix[each_atom][1:len(adj_idx_matrix)]])

            if adj_ele_matrix[each_atom][0] == 'O':
                ord_adj_ele_matrix.append(['O',adj_ele_matrix[each_atom][1:len(adj_ele_matrix)]])
                ord_adj_idx_matrix.append([each_atom+1,adj_idx_matrix[each_atom][1:len(adj_idx_matrix)]])

            if adj_ele_matrix[each_atom][0] == 'N':
                ord_adj_ele_matrix.append(['N',adj_ele_matrix[each_atom][1:len(adj_ele_matrix)]])
                ord_adj_idx_matrix.append([each_atom+1,adj_idx_matrix[each_atom][1:len(adj_idx_matrix)]])

            if adj_ele_matrix[each_atom][0] == 'C':
                ord_adj_ele_matrix.append(['C',adj_ele_matrix[each_atom][1:len(adj_ele_matrix)]])
                ord_adj_idx_matrix.append([each_atom+1,adj_idx_matrix[each_atom][1:len(adj_idx_matrix)]])

            if adj_ele_matrix[each_atom][0] == 'H':
                ord_adj_ele_matrix.append(['H',adj_ele_matrix[each_atom][1:len(adj_ele_matrix)]])
                ord_adj_idx_matrix.append([each_atom+1,adj_idx_matrix[each_atom][1:len(adj_idx_matrix)]])

        #now order all the atoms in each branch
        for each_atom in range(n_targets):
        
            neighbors = ord_adj_ele_matrix[each_atom][1]

            ord_nbridx_list = []
            ord_nbrele_list = []


            #get fluorines first
            for each_nbr in range(len(neighbors)):
                if neighbors[each_nbr] == 'F':
                    ord_nbrele_list.append('F')
                    ord_nbridx_list.append(ord_adj_idx_matrix[each_atom][1][each_nbr])
            #get oxygens second
            for each_nbr in range(len(neighbors)):
                if neighbors[each_nbr] == 'O':
                    ord_nbrele_list.append('O')
                    ord_nbridx_list.append(ord_adj_idx_matrix[each_atom][1][each_nbr])

            for each_nbr in range(len(neighbors)):
                if neighbors[each_nbr] == 'N':
                    ord_nbrele_list.append('N')
                    ord_nbridx_list.append(ord_adj_idx_matrix[each_atom][1][each_nbr])

            for each_nbr in range(len(neighbors)):
                if neighbors[each_nbr] == 'C':
                    ord_nbrele_list.append('C')
                    ord_nbridx_list.append(ord_adj_idx_matrix[each_atom][1][each_nbr])

            for each_nbr in range(len(neighbors)):
                if neighbors[each_nbr] == 'H':
                    ord_nbrele_list.append('H')
                    ord_nbridx_list.append(ord_adj_idx_matrix[each_atom][1][each_nbr])

            #replace ordered for the one in the original list
            ord_adj_ele_matrix[each_atom][1] = ord_nbrele_list
            ord_adj_idx_matrix[each_atom][1] = ord_nbridx_list
    
    
    if n_targets == 1:

        ord_adj_idx_matrix = adj_idx_matrix
        ord_adj_ele_matrix = adj_ele_matrix

        neighbors = ord_adj_ele_matrix

        ord_nbridx_list = []
        ord_nbrele_list = []

        #get fluorines first
        for each_nbr in range(len(neighbors)):
            if neighbors[each_nbr] == 'F':
                ord_nbrele_list.append('F')
                ord_nbridx_list.append(ord_adj_idx_matrix[each_nbr])
        #get oxygens second
        for each_nbr in range(len(neighbors)):
            if neighbors[each_nbr] == 'O':
                ord_nbrele_list.append('O')
                ord_nbridx_list.append(ord_adj_idx_matrix[each_nbr])

        for each_nbr in range(len(neighbors)):
            if neighbors[each_nbr] == 'N':
                ord_nbrele_list.append('N')
                ord_nbridx_list.append(ord_adj_idx_matrix[each_nbr])

        for each_nbr in range(len(neighbors)):
            if neighbors[each_nbr] == 'C':
                ord_nbrele_list.append('C')
                ord_nbridx_list.append(ord_adj_idx_matrix[each_nbr])

        for each_nbr in range(len(neighbors)):
            if neighbors[each_nbr] == 'H':
                ord_nbrele_list.append('H')
                ord_nbridx_list.append(ord_adj_idx_matrix[each_nbr])
        
        ord_adj_ele_matrix = ord_nbrele_list
        ord_adj_idx_matrix = ord_nbridx_list

    return ord_adj_idx_matrix, ord_adj_ele_matrix

def depth_label(adj_idx_matrix,adj_ele_matrix,target_connect_list,target_element_list,n_depth):
    '''
    Main chunk of the code that runs through targets/roots of a molecule 
    to update the connection_matrix, element_matrix for this depth

        n_branch                - initialize the number of branches for each atom/row of adjacency 
                                  which is just automatically 1, this will change as atom picks up connections and branches
        target_connect_list     - this is modified to be a list of list ([target_connect_list]), to allow for branching at the target/root/atom. 
                                  Initially at the root, there will be no branching, but as the depths pick up, there might be branching 
                                  in which we include a list for each branch
                                    
                                    [[ root, [root_connections]]_BRANCH]
        
        save_root               - in order not to go backwards through the depths, we save all roots/atoms as a list we have already done
                                  and make sure they are never included again as neighbors. N
        
        target_depth_idxs       - this will hold the list of idxs of all the branches at all the depths, it adds another bracket/list for the depths
                                  
                                    [ (depth) [ (branch) [root, [root_connections]] ] ]

        target_depth_elements   -this will hold the list of elemends of all the branches at all the depths, it adds another bracket/list for the depths
                                  
                                    [ (depth) [ (branch) [root, [root_connections]] ] ]
                                    
        n_depth                 - the specified depth for the label 
        nexttarget_list         - this is at first emppty but as we go through the neighbors of each branch of a previous depth, 
                                  it finds all the branches of the next depth and also gives information about the depth after that (the neighbors of the roots of the branches)
                                  all the branches found from the new depth's branches are added to target_depth_idxs 
                                  as well to keep track of all the depths
        nextelement_list        - this is also at first empty, but as we go through the neighbors of a previous depth's branch,
                                  will hold the branch elements at this depth, and has information about next branch elements (same as above but for elements)
        each_branch             - an integer that that runs through all the branches. initially there is just one branch at the root atom. 
                                  but as the depth picks up, there may be multiple of neighbors that will each get a list of their own as branches 
                                  as defined above and ALL THE BRANCHES FOR A DEPTH WILL BE HELD in new_target_list before being cleared for the next depth
        nbr_list                - the list of neighbors for each branch's root/atom
        n_neighbor              - number of neighbors in the branch
        each_neighbor           - an integer that runs through all the neighbors in the branch
        nbr                     - singling out a neighbor from the branch, 
                                  saving that neighbor in save_roots
                                  find the neighbors' neighbors using the adjacency matrix 
                                  and then append that as a new branch in next_target_list
        each_connection         - an integer that runs through the adjacency matrix just to find the nbr and
                                  get the neighbors of that neighbor as the new target branches. The nbr now becomes the root
                                  together with the nbrs of nbr will be another branch added to the next_target_list and another branch added to the full_connect_list as well
                                  EVERY nbr will add their own branch of nbrs-nbrs UNLESS we are going backwards NO going backwards
        target_depth_idxs        - appends all the branches of idxs (nbrs) at each depth which are found in next_target_list
        target_depth_elements    - appends all the branches of elements (nbrs) at each depth which are found in next_target_list
    '''    
    
    #Initialize # of branches for each target passing
    #Which is always 1 for each target passing,
    #this will change as we run depth labelling on each target separately
    #[[]] means 1 branch, [[],[]] means two branches of a bunch of connections (thus new targets) each
    n_branch = 1
    #store target_connect_list as having a single branch to start
    target_connect_list = [target_connect_list]  #NOW we have [[ root, [root_connections]]_branch]
    target_element_list = [target_element_list]

    #In order not to go backwards and keep pushing forward we should keep track
    #of all roots we have already passed and make sure not to pass them again
    #NOTE DO NOT MESS AROUND WITH THE BRACKETS
    save_roots = [target_connect_list[0][0]]
        
    #Define a full connection matrix which will hold all the connecting indicies
    #THROUGH DEPTH of target
    #NOTE that this has to be defined as A LIST OF LISTS
    #This implies there is a starting branch of 1
    target_depth_idxs = [target_connect_list] # now we have [ [ [ (root) , [ root_connections ] ]_BRANCH]_DEPTH ] 
    target_depth_elements = [target_element_list]


    #if n_depth is zero, then this is just the roots of each of the target branch passing through
    #that is the atoms of the molecule, 
    #LATER ON THIS IS HANDLED DURING THE FG_KEY 
    #if n_depth == 1, then this has already been done because
    #adjacency already does 1 depth and this can be reduced to simply
    #returning target_depth_elements and working off of that 1 depth - 1 branch
   
    if n_depth == 0 or n_depth == 1:
        target_depth_elements = target_depth_elements
    #If total specified depth > 1, we go through each target's connections
    #And start branching from that target
    elif n_depth > 1:
        # run through the rest of the depth (except 1 because we already did the first one)
        for each_depth in range(n_depth-1):

            #define a nexttarget_list which will hold the next targets
            nexttarget_list = []
            nextelement_list = []

            for each_branch in range(n_branch):

                #define neighbor list
                nbr_list = target_connect_list[each_branch][1]
                n_neighbor = len(nbr_list)

                for each_neighbor in range(n_neighbor):
                    
                    nbr = target_connect_list[each_branch][1][each_neighbor]

                    if nbr not in save_roots:
                        save_roots.append(nbr)

                        #if neighbor is not in saved roots
                        #find each neighbor's connect_list and element_list and store them in a set of lists
                        for each_connection in range(len(adj_idx_matrix)):
                            if adj_idx_matrix[each_connection][0] == nbr:
                                nexttarget_list.append(adj_idx_matrix[each_connection])
                                nextelement_list.append(adj_ele_matrix[each_connection])

            #append neighbor's connections to new_connect list
            target_depth_idxs.append(nexttarget_list)
            target_depth_elements.append(nextelement_list)

            #newconnect_list now becomes the "old" conncet_list and we try again for the next n_depth
            target_connect_list = nexttarget_list
            target_element_list = nextelement_list

            n_branch = len(target_connect_list)

    return target_depth_elements, target_depth_idxs


def order_branches_in_each_depth(target_depth_elements,target_depth_connections,n_depth):
    '''
    This is the second round of ordering the target_depth_elements
    It is meant to ensure that the branches are odered in a consistent way
    (BECAUSE THERE MAY BE (AND OFTEN IS) MULTIPLE IN EACH DEPTH THAT HAVE SAME ROOT ATOM)
    
    This is done using a rating, which adds all
    the atomic number priority element-rating in a branch
    and compares the value of the rating to other branches in the depth
    and then orders the branches in each depth accordingly
    This ensures that random order of listing branches in the fg_key 
    is the SAME and consistent throughout

        each_depth               - an integer that runs through all the depths of the target branches collected
        target_depth_elements    - the result of the depth_label on the ordered adjacency/branches elements
        target_depth_connections - the results of the depth_label on the ordered adacency/branches atom idxs
        n_branch                 - number of branches for each depth
        brch                     - branches of each depth
        rating                   - sums the rating of each branch based on atomic number priority
                                   'F' gets exponentially higher rating than 'O' always ensuring that 
                                   if one 'F' is found then that branch is first, if 2 then even better... etc
        newbranch_order          - empty list initially that will hold all branches in the correct order 
                                   according to the rating given to each
        each_branch              - integer running through all the branches at a depth

        n_neighbor               - number of neighbors in a branch
        each_neighbor            - integer running throuhg all the neighbors in a branch
        score                    - sums the value of all neighbors in a branch 
                                   to append the score of that branch to the ratings list of the branches
        rating_copy              - keep a copy of rating list because we will be using pop on the copy,
                                   and perserving the original, we will remove the max_value index each time on the copy
                                   adding it to the ordered branch list (newbranch_order) 
                                   and repeating until we are done with all the branches 
        max_value                - finds the branch with the maximum rating 
        max_index                - finds the index of the branch with the maximum rating 
                                   (if two branches have equal rating then it is okay if they switch around because that means they HAVE equivalent neighbors)
                                   uses that index to append the correct branch index to the newbranch_order list
        max_index_copy           - if you use just the max_index to pop rating_copy, you will run into the problem of 
                                   max_index of rating not matching the max_index of rating_copy because you deleted some of rating_copy
                                   so when you find max_index of the rating (using max_value of rating copy) if will be a DIFFERENT max_index
                                   than the max_index_copy FOR rating_copy
    '''

    

    for each_depth in range(len(target_depth_elements)):

        #define #number of branches
        n_branch = len(target_depth_elements[each_depth])
        brch_ele = target_depth_elements[each_depth]
        brch_con = target_depth_connections[each_depth]

        #if there are more than one branch, 
        #give a rating to each branch and reorder the branches according to rating
        rating = []
        newbranchele_order = []
        newbranchcon_order = []
        if n_branch == 1:
            newbranchele_order = target_depth_elements[each_depth]
            newbranchcon_order = target_depth_connections[each_depth]

        if n_branch > 1:
            
            for each_branch in range(n_branch):
            #for each neighbor list sum priority numbers
                n_neighbor = len(target_depth_elements[each_depth][each_branch][1])

                score = 0
                #add score of the root atom
                if target_depth_elements[each_depth][each_branch][0] == 'F':
                    score = score + 1000000 + 1/(each_branch+1)
                if target_depth_elements[each_depth][each_branch][0] == 'O':
                    score = score + 100000 + 1/(each_branch+1)
                if target_depth_elements[each_depth][each_branch][0] == 'N':
                    score = score + 10000 + 1/(each_branch+1)
                if target_depth_elements[each_depth][each_branch][0] == 'C':
                    score = score + 1000 + 1/(each_branch+1)
                if target_depth_elements[each_depth][each_branch][0] == 'H':
                    score = score + 100 + 1/(each_branch+1)


                for each_neighbor in range(n_neighbor):

                    if target_depth_elements[each_depth][each_branch][1][each_neighbor] == 'F':
                        score = score + 1000000 + 1/(each_branch+1)
                    if target_depth_elements[each_depth][each_branch][1][each_neighbor] == 'O':
                        score = score + 100000 + 1/(each_branch+1)
                    if target_depth_elements[each_depth][each_branch][1][each_neighbor] == 'N':
                        score = score + 10000 + 1/(each_branch+1)
                    if target_depth_elements[each_depth][each_branch][1][each_neighbor] == 'C':
                        score = score + 1000 + 1/(each_branch+1)
                    if target_depth_elements[each_depth][each_branch][1][each_neighbor] == 'H':
                        score = score + 100 + 1/(each_branch+1)
                rating.append(score)
            
                #go through max rating index as index for ordering each branch
            index_done = []
            rating_copy = rating.copy()
            for each_branch in range(n_branch):
                max_value = max(rating_copy)

                max_index = rating.index(max_value)

                newbranchele_order.append(brch_ele[max_index])
                newbranchcon_order.append(brch_con[max_index])
                
                max_index_copy = rating_copy.index(max_value)
                rating_copy.pop(max_index_copy)

        #replace n_depth
        target_depth_elements[each_depth] = newbranchele_order
        target_depth_connections[each_depth] = newbranchcon_order

    return target_depth_elements, target_depth_connections

def fg_key(target_depth_elements,n_depth):
    '''
    The final result of the entire process is making the fg_key for the atomic environment found
    through the depths (the arrays target_depth_idxs, target_depth_elements)
    this is directly based off of our target_depth_element list, 
    we just replace the brackets in each depth/branch/neighborlist to make the fg_key


        n_depth         - the depth of the label specified by the user
                             if the depth == 0 then this is just the target/root atoms
                             NOT the branches

                             else --> we build the fg_key for the branches at each depth
        
        fg_key          - the string label that will identify this distinct depth atomic environment from the rest
        
                            it follows this format

                            'rootelement-[branching_nbrs]-rootelement-[branching_nbrs]'
                            
                            where rootelement are the roots of the branches
                            and the branching_nbrs are the "leaves" or "tips" of the branches

                            the above will always be ordered according to ATOMIC NUMBER PRIORTY

        each_depth      - an integer that runs through all the depths of target_depth_elements
        n_branch        - number of branchs at each depth
        brch            - all the branches at each depth
        each_branch     - integer running through the number of branches at a depth. Places
                          the root brch[each_branch][0] as 'root-'
        n_neighbor      - number of neighbors in each branch
        each_neighbor   - integer running through index of each neighbor up to 
                          number of neighbors in a branch. Places all neighbors in
                          '[neighbors]-' in the fg_key. 
    '''

    if n_depth == 0:
        fg_key = target_depth_elements[0][0][0]
    #just order the branches according to most oxygen,nitrogen,hydrogen
    else:
        #construct fg key directly from organized list
        fg_key = ''
        for each_depth in range(len(target_depth_elements)):
#            fg_key = fg_key + '-'

            n_branch = len(target_depth_elements[each_depth])
            brch = target_depth_elements[each_depth]

            for each_branch in range(n_branch):

                fg_key = fg_key + brch[each_branch][0] + '-'

                n_neighbor = len(target_depth_elements[each_depth][each_branch][1])
                
                fg_key = fg_key + '['

                for each_neighbor in range(n_neighbor):

                    fg_key = fg_key + target_depth_elements[each_depth][each_branch][1][each_neighbor] + '-'
                fg_key = fg_key + ']-'
            
            fg_key = fg_key + '|-'
    return fg_key


def get_emb_avgs(embeddings,n_features):
    '''
    calculates average embeddings per atom type and returns as a dictionary to each atom-type

    inputs:
        embeddings          - embeddings of the entire data labelled according to atom_num
        n_features          - number of features per atom-embedding

    processing:
        embs_Z              - dataframe boolean indexing to find all embeddings with atom_num = Z
        embs_Zavg           - taking mean of the dataframe for the embeddings of atom Z

    outputs:
        embs_avg_dict       - the dictionary that matches atom-type to embedding average for that atom-type
                              returned to be used to build the zeroth order perturbation in terms of embeddings
    
    '''

    embs_H = embeddings[embeddings['atom_num'] == 1].values
    embs_Havg = np.mean(embs_H[:,0:n_features],axis=0)

    embs_C = embeddings[embeddings['atom_num'] == 6].values
    embs_Cavg = np.mean(embs_C[:,0:n_features],axis=0)

    embs_N = embeddings[embeddings['atom_num'] == 7].values
    embs_Navg = np.mean(embs_N[:,0:n_features],axis=0)

    embs_O = embeddings[embeddings['atom_num'] == 8].values
    embs_Oavg = np.mean(embs_O[:,0:n_features],axis=0)

    embs_F = embeddings[embeddings['atom_num'] == 9].values
    embs_Favg = np.mean(embs_F[:,0:n_features],axis=0)

    #make a dictionary of embs5 averages for each element so that they are easy to get
    embs_avg_dict = {'H': embs_Havg,
                    'C': embs_Cavg,
                    'N': embs_Navg,
                    'O': embs_Oavg,
                    'F': embs_Favg}
    
    return embs_avg_dict


def pertembs(n_depth,target_depth_connections, target_depth_elements, embedings_mol,embs_avg_dict, n_features, n_pert, X, Y, target_idx, prev_embs_mol):
    '''
    sets up the linear perturbation problem using the idx label
    
    
    '''


    depth_embs = []
    for depth_i, depth_atomidxs in enumerate(target_depth_connections):

        number_max_branches_atdepth = 4**(depth_i)

        branch_embs = []
        count_branches = 0
        for branch_i, branch in enumerate(target_depth_connections[depth_i]):
            
            branch_elements = target_depth_elements[depth_i][branch_i]

            neighbor_embs = []
            count_nbrs = 0
            for neighbor_i, neighbor in enumerate(branch[1]):

                neighbor_ele = branch_elements[1][neighbor_i]
                
                if n_pert == 1:
                    neighbor_emb = embs_avg_dict[neighbor_ele]
                elif n_pert > 1:
                    neighbor_emb = prev_embs_mol[neighbor-1] #-1 to match the idx of python because these idxs come from mol files
                else:
                    raise ValueError('n_pert cannot be less than 1')

                neighbor_embs = np.hstack((neighbor_embs,neighbor_emb))

                count_nbrs = count_nbrs + 1

            if count_nbrs < 4:
                missing_nbrs = 4 - count_nbrs
                neighbor_embs = np.hstack((neighbor_embs,np.zeros((n_features*missing_nbrs))))

            branch_embs = np.hstack((branch_embs,neighbor_embs))
            count_branches =  count_branches + 1

        if count_branches < number_max_branches_atdepth:
            missing_branches = number_max_branches_atdepth - count_branches
            branch_embs = np.hstack((branch_embs,np.zeros(missing_branches*n_features*4)))
    
        depth_embs = np.hstack((depth_embs,branch_embs))

    Y = np.vstack((Y,embedings_mol.iloc[target_idx]))
    X = np.vstack((X,depth_embs))


    return X, Y