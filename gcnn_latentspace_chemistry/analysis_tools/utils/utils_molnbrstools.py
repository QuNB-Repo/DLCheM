#A function that helps you get the id of an atom given it's atomindex
def get_atom_id(atomindex,mol_filepath):
    mol_string_filewrapper = open(mol_filepath,mode='r')
    for line_index, each_line in enumerate(mol_string_filewrapper):
        if line_index == 3 + atomindex:
            nbr_id = str(each_line[31:32])

    return nbr_id

#get indices of each connected atom, i.e neighbors 
def find_nextdoor_neighbors(depth,mol_filepath,atomidx,number_atoms):
    
    #open the molecule's mol file
    mol_string_filewrapper = open(mol_filepath,mode='r')
    
    #initialize an empty list that will store all nbr indexs and nbr identities
    nbrs = []
    nbrs_ids = []
    
    #if depth is only one
    if depth == 1: 

        #Go through each line of the mol file
        for line_index, each_line in enumerate(mol_string_filewrapper): 

            #After 4 lines + number atoms comes the bonding section
            #This will help us find the neighbors (indices and ids) bonded directly
            if line_index > 3 + number_atoms:

                #look for target idx in the first 7 characters of each line only, for no conflicts, (this may be redundant line)
                if str(atomidx) in each_line[0:7]:

                    #if it is in the first 4 characters, then the nbr is in the next 4
                    #NOTE that the target index as string should be looked for as ' i ' with spaces on either side
                    #so that double digits numbers are not also included in the count, i.e 'ii'
                    if ' ' +str(atomidx)+' ' in each_line[0:4]:
                        nbr = int(each_line[4:7])
                        nbr_id = get_atom_id(nbr,mol_filepath)
                        nbrs.append(nbr)
                        nbrs_ids.append(nbr_id)

                    #if it is in the second 4 characters, then the nbr is in the first 4
                    #same NOTE as above
                    if ' ' +str(atomidx)+' ' in each_line[3:7]:
                        nbr = int(each_line[0:4])
                        nbrs.append(nbr)
                        nbr_id = get_atom_id(nbr,mol_filepath)
                        nbrs_ids.append(nbr_id)

    #NOT yet implemented for greater depths, need to talk to stijn
    else:
        return NotImplementedError
    

    #get the neighbors identity 
    return nbrs, nbrs_ids