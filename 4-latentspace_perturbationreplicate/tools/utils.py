import numpy as np
import fileinput
import os

# Set a seed for reproducibility
np.random.seed(42)

def find_nextdoor_neighbors(MOL_FILEPATH,ATOMIDX,NUMBER_ATOMS):
        ''' 
        looks up an atom's direct neighbors using the mol file and the atom's index

            Args:
                MOL_FILEPATH        ---- where the temporary mole file is located for the molecule
                ATOMIDX             ---- which atom in the molecule, by index, is the target atom for finding next door neighbors
                NUMBER_ATOMS        ---- number of atoms in the molecule for help with reading the mol file

        '''
        #mol_string_file ---- opening mol file in read mode
        mol_string_file = open(MOL_FILEPATH,mode='r')

        #A quick function that helps you get the id of an atom (NEIGHBOR) given it's atomindex
        def get_atom_id(nbridx,mol_filepath):
            #mol_string_filewrapper ---- the mol filepath has to be opened again in read mode, I found this not to work if I just just use the mol string in here, still cannot figure out why (minor)
            mol_string_filewrapper = open(mol_filepath,mode='r')
            #for each line in the mol string
            for line_index, each_line in enumerate(mol_string_filewrapper):
                #if line is at the location of the neighbor's index
                if line_index == 3 + nbridx:
                    #nbr_id ---- neighbor's element is found on line 31.
                    nbr_id = str(each_line[31:32])

            return nbr_id
        
        #initialize an empty list that will store all nbr indexs and nbr identities
        #nbrs ---- will store the neighbor idxs of the target 
        #nbrs ---- will store the neighbor elements of the target
        #bond_nmbrs ---- will store the bond order with the neighbors of target
        nbrs = []
        nbrs_ids = []
        bond_nmbrs = []

        #Go through each line of the mol file, enumerate it as well so we count which line index we are on
        for line_index, each_line in enumerate(mol_string_file): 

            #After 4 lines + number atoms comes the bonding section
            #This will help us find the neighbors (indices and ids) bonded directly
            if line_index > 3 + NUMBER_ATOMS:

                #look for target idx in the first 7 characters of each line only, for no conflicts, (this may be redundant line)
                if str(ATOMIDX) in each_line[0:7]:

                    #if it is in the first 4 characters, then the nbr is in the next 4
                    #NOTE that the target index as string should be looked for as ' i ' with spaces on either side
                    #so that double digits numbers are not also included in the count, i.e 'ii'
                    if ' ' +str(ATOMIDX)+' ' in each_line[0:4]:
                        #nbr ---- extracting neighbor as the 4th to 7th character in the line
                        nbr = int(each_line[4:7])
                        #nbr_id ---- getting the neighbor element using our get_atom_id function now that we know nbr's index
                        nbr_id = get_atom_id(nbr,MOL_FILEPATH)
                        #nbrs ---- append nbrs to nbrs list
                        nbrs.append(nbr)
                        #nbrs_ids ---- appends nbr element to nbrs_ids list
                        nbrs_ids.append(nbr_id)

                        #number_bonds ---- get number of bonds for the atom, found in characters 7:9
                        number_bonds = int(each_line[7:9])
                        #bond_nmbrs ---- append nbr bond order to bond_nmbrs lists
                        bond_nmbrs.append(number_bonds)

                    #if it is in the second 4 characters, then the nbr is in the first 4
                    #NOTE that the target index as string should be looked for as ' i ' with spaces on either side
                    #so that double digits numbers are not also included in the count, i.e 'ii'
                    if ' ' +str(ATOMIDX)+' ' in each_line[3:7]:
                        #nbr ---- extracting neighbor as the 0th to 4th character in the line
                        nbr = int(each_line[0:4])
                        #nbrs ---- append nbrs to nbrs list
                        nbrs.append(nbr)
                        #nbr_id ---- getting the neighbor element using our get_atom_id function now that we know nbr's index
                        nbr_id = get_atom_id(nbr,MOL_FILEPATH)
                        #nbrs_ids ---- appends nbr element to nbrs_ids list
                        nbrs_ids.append(nbr_id)
                
                        #number_bonds ---- get number of bonds for the atom, found in characters 7:9
                        number_bonds = int(each_line[7:9])
                        #bond_nmbrs ---- append nbr bond order to bond_nmbrs lists
                        bond_nmbrs.append(number_bonds)
                        
        #get the neighbors identity 
        return nbrs, nbrs_ids, bond_nmbrs


def convertxyz_to_molfull(ATOMICNUMBERS,ATOMICPOSITIONS,TEMP_FILENAME):
    '''
    function that constructs xyz file and convert to mol file using atomic number and atomic positions information
    
        Args:
            ATOMICNUMBERS       - numpy row array of the atomic numbers in the molecule 
            ATOMICPOSITIONS     - numpy matrix array of the x,y,z positions of each atom 
            TEMP_FILENAME       - temporary filename to save the xyz and convert to mol file with obabel

        Returns:
            tempmolfilepath     - mol file name that was constructed for use by algorithm in finding neighborhood atoms    
    '''

    #atomic_num_to_symbols ---- dictionary converting atomic number to dictionary for writing xyz
    atomic_num_to_symbols = {
            1: "H",
            6: "C",
            7: "N",
            8: "O",
            9: "F"
        }

    #element_symbols ---- converting atomic numbers array to atomic symbols array using dictionary
    element_symbols = [atomic_num_to_symbols[num] for num in ATOMICNUMBERS]
    
    #num_atoms ---- integer number of atoms in the molecule
    num_atoms = len(ATOMICNUMBERS)
    
    #xyz_string ---- wrting the xyz string for then writing onto an xyz file
    xyz_string = f"{num_atoms}\n0.00000\n"

    #for each element_symbol, and x,y,z, in element_symbols, ATOMICPOSITIONS, use to write each line in xyz_string 
    for element_symbol, (x, y, z) in zip(element_symbols, ATOMICPOSITIONS):
        #xyz_string ---- updated with each line of element sumbol and coordinates
        xyz_string += f"{element_symbol} {x:.5f} {y:.5f} {z:.5f}\n"

    #tempxyzfilepath ---- temporary xyz filepath to write xyz to so that we can convert to a mol file
    tempxyzfilepath = TEMP_FILENAME+".xyz"
    #open tempxyzfilepath and write the xyz string
    with open(tempxyzfilepath, 'w') as file:
        #file ---- writing xyz string onto temp file
        file.write(xyz_string)

    #tempmolfilepath ---- temporary mol filepath to convert from xyz, will be placed in the same directory 
    tempmolfilepath = TEMP_FILENAME+".mol"
    #out ---- output of an os.system operation which runs obabel from command line to convert xyz to mol file
    out = os.system('obabel ' + tempxyzfilepath +' -O ' + tempmolfilepath + '')

    #return the name of the mol filepath so it can be found later by the algorithm
    return tempmolfilepath



def fill_nbrhood_featurespace_minimal(NBRHOOD_EMBED_FEATURESPACE,NBRS_IDXS,NBRS_IDS,DEPTH_i,N_FEATURES,PREV_Y_MOL):
    '''
    fills the nbrhood_embed_featurespace for each atom with the embedding of its neighbors from previous perturbation, 
    in the proper element-embedding and neighborhood-depth place holder 

        Args:
            NBRHOOD_EMBED_FEATURESPACE      - each atom's feature space that is used to predict peturbational corrections,
                                              it is made up of 5*DEPTH*EMB_FEATURES features , where 5 represents element-embedding (EMB_FEATURES) placeholds for H,C,N,O,F for each DEPTH
            NBRS_IDXS                       - the list of nbr indices for this atom
            NBRS_IDS                        - the list of nbr elements for this atom
            DEPTH_i                         - which depth of neighborhood to fill (i.e which layer of neighborhood for this atom is extracted at current state) 
            N_FEATURES                      - number of features per atom-embedding
            PREV_Y_MOL                      - prev perturbational correction results which will now be used to fill the nbrhood feature space
        Returns:
            NBRHOOD_EMBED_FEATURESPACE      - updated so that this depth of neighbors are filled
    '''
    #for each neighhbor for this atom
    for nbr_i in range(len(NBRS_IDXS)):

        #nbr_idx ---- converting from mol indexing (starts at 1) to python indexing (starts at 0) for each neighbor index, 
        #used to pick out prev Y perturbation for the neighbor atom
        nbr_idx = NBRS_IDXS[nbr_i]-1

        #if H neighbor
        if NBRS_IDS[nbr_i] == 'H': 
            #fill H-embedding placeholder (1) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[DEPTH_i*5*N_FEATURES:(DEPTH_i*5+1)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[DEPTH_i*5*N_FEATURES:(DEPTH_i*5+1)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]
        #if C neighbor
        if NBRS_IDS[nbr_i] == 'C':
            #fill C-embedding placeholder (2) with prev_y  at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+1)*N_FEATURES:(DEPTH_i*5+2)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+1)*N_FEATURES:(DEPTH_i*5+2)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]        
        #if N neighbor       
        if NBRS_IDS[nbr_i] == 'N': 
            #fill N-embedding placeholder (3) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+2)*N_FEATURES:(DEPTH_i*5+3)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+2)*N_FEATURES:(DEPTH_i*5+3)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]
        #if O neighbor   
        if NBRS_IDS[nbr_i] == 'O':
            #fill O-embedding placeholder (4) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+3)*N_FEATURES:(DEPTH_i*5+4)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+3)*N_FEATURES:(DEPTH_i*5+4)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]
        #if F neighbor  
        if NBRS_IDS[nbr_i] == 'F': 
            #fill F-embedding placeholder (5) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+4)*N_FEATURES:(DEPTH_i*5+5)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*5+4)*N_FEATURES:(DEPTH_i*5+5)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]   

    #return NBRHOOD_EMBED_FEATURESPACE for this atom after all its neighbors embedding from prev perturbation have been filled
    return NBRHOOD_EMBED_FEATURESPACE


def fill_nbrhood_featurespace_extended(NBRHOOD_EMBED_FEATURESPACE,NBRS_IDXS,NBRS_IDS,NBRS_BOND_NMBRS,DEPTH_i,N_FEATURES,PREV_Y_MOL):
    '''
    fills the nbrhood_embed_featurespace for each atom with the embedding of its neighbors from previous perturbation, 
    in the proper element-embedding and neighborhood-depth place holder

    this time its been extended to include bond order with neighbor as part of the feature space 

        Args:
            NBRHOOD_EMBED_FEATURESPACE      - each atom's feature space that is used to predict peturbational corrections,
                                              it is made up of 5*DEPTH*EMB_FEATURES features , where 5 represents element-embedding (EMB_FEATURES) placeholds for H,C,N,O,F for each DEPTH
            NBRS_IDXS                       - the list of nbr indices for this atom
            NBRS_IDS                        - the list of nbr elements for this atom
            NBRS_BOND_NMBRS                 - the list of nbr bond orders for this atom
            DEPTH_i                         - which depth of neighborhood to fill (i.e which layer of neighborhood for this atom is extracted at current state) 
            N_FEATURES                      - number of features per atom-embedding
            PREV_Y_MOL                      - prev perturbational correction results which will now be used to fill the nbrhood feature space
        Returns:
            NBRHOOD_EMBED_FEATURESPACE      - updated so that this depth of neighbors are filled
    '''

    #for each neighhbor for this atom
    for each_nbr in range(len(NBRS_IDXS)):

        #nbr_idx ---- converting from mol indexing (starts at 1) to python indexing (starts at 0) for each neighbor index, 
        #used to pick out prev Y perturbation for the neighbor atom
        nbr_idx = NBRS_IDXS[each_nbr]-1

        #if H neighbor
        if NBRS_IDS[each_nbr] == 'H': 
            #fill H-embedding placeholder (1) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[DEPTH_i*10*N_FEATURES:(DEPTH_i*10+1)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[DEPTH_i*10*N_FEATURES:(DEPTH_i*10+1)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]
        #if C neighbor with bond order 1    
        if NBRS_IDS[each_nbr] == 'C' and NBRS_BOND_NMBRS[each_nbr] == 1: 
            #fill C-embedding bond order 1 placeholder (2) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+1)*N_FEATURES:(DEPTH_i*10+2)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+1)*N_FEATURES:(DEPTH_i*10+2)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]        
        #if C neighbor with bond order 2        
        if NBRS_IDS[each_nbr] == 'C' and NBRS_BOND_NMBRS[each_nbr] == 2: 
            #fill C-embedding bond order 2 placeholder (3) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+2)*N_FEATURES:(DEPTH_i*10+3)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+2)*N_FEATURES:(DEPTH_i*10+3)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]                                                               
        #if C neighbor with bond order 3
        if NBRS_IDS[each_nbr] == 'C' and NBRS_BOND_NMBRS[each_nbr] == 3: 
            #fill C-embedding bond order 2 placeholder (4) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+3)*N_FEATURES:(DEPTH_i*10+4)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+3)*N_FEATURES:(DEPTH_i*10+4)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]  

        #if N neighbor with bond order 1        
        if NBRS_IDS[each_nbr] == 'N' and NBRS_BOND_NMBRS[each_nbr] == 1: 
            #fill N-embedding bond order 1 placeholder (5) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+4)*N_FEATURES:(DEPTH_i*10+5)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+4)*N_FEATURES:(DEPTH_i*10+5)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]
        #if N neighbor with bond order 2
        if NBRS_IDS[each_nbr] == 'N' and NBRS_BOND_NMBRS[each_nbr] == 2: 
            #fill N-embedding bond order 2 placeholder (6) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+5)*N_FEATURES:(DEPTH_i*10+6)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+5)*N_FEATURES:(DEPTH_i*10+6)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]
        #if N neighbor with bond order 3
        if NBRS_IDS[each_nbr] == 'N' and NBRS_BOND_NMBRS[each_nbr] == 3: 
            #fill N-embedding bond order 3 placeholder (7) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+6)*N_FEATURES:(DEPTH_i*10+7)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+6)*N_FEATURES:(DEPTH_i*10+7)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]                        

        #if O neighbor with bond order 1
        if NBRS_IDS[each_nbr] == 'O' and NBRS_BOND_NMBRS[each_nbr] == 1:
            #fill O-embedding bond order 1 placeholder (8) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+7)*N_FEATURES:(DEPTH_i*10+8)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+7)*N_FEATURES:(DEPTH_i*10+8)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]
        #if O neighbor with bond order 2       
        if NBRS_IDS[each_nbr] == 'O' and NBRS_BOND_NMBRS[each_nbr] == 2: 
            #fill O-embedding bond order 2 placeholder (9) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+8)*N_FEATURES:(DEPTH_i*10+9)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+8)*N_FEATURES:(DEPTH_i*10+9)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]                           

        #if F neighbor
        if NBRS_IDS[each_nbr] == 'F': 
            #fill F-embedding bond order 1 placeholder (10) with prev_y at the specified depth
            NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+9)*N_FEATURES:(DEPTH_i*10+10)*N_FEATURES] = NBRHOOD_EMBED_FEATURESPACE[(DEPTH_i*10+9)*N_FEATURES:(DEPTH_i*10+10)*N_FEATURES] + PREV_Y_MOL[nbr_idx][0:N_FEATURES]   

    #return NBRHOOD_EMBED_FEATURESPACE for this atom after all its neighbors embedding from prev perturbation have been filled
    return NBRHOOD_EMBED_FEATURESPACE




