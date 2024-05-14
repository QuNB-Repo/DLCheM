from ase.io import read
import numpy as np
from schnetpack import AtomsData

import os

import warnings 

warnings.filterwarnings("ignore")

'''
Last Updated: 2024-03-06


The various tools used in FGTransform and analysis

    generate_dbfromxyz  - generates db files from xyz so they can be accessed by SchNet models
                          (for embedding extractions)
    convert_xyz2mol     - converts xyz files to mol files (important for functional group bonding)

'''

def generate_dbfromxyz(XYZDATASET_FILEPATH,AVAILABLE_PROPERTIES):
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
            XYZDATASET_FILEPATH         - the xyz dataset filepath to be converted to db
            AVAILABLE_PROPERTIES        - list of strings representing the names of the properties

        Process:
            number_properties           - number of properties for each molecule
            molecules                   - uses the atomic simulation environment (ase) package to read 
                                          the xyz dataset filepath and ensure that it is separated and indexed
                                          by molecular 'frames'
            properties_list             - list that will hold all the properties for ALL molecules for the db conversion
                                          each molecule will get its own dictionary of properties
                                          with the property name being the key and the value being the property itself

                                           [{prop1: #value1, prop1: #value2}_mol1,{}_mol2... {}_molN] for each molecule

            properties                  - a list that holds all the properties for each molecule to be added to properties_list
            mol                         - each molecule in the xyz ase read file
            each_prop                   - index that runs over all the properties available
            db_dataset_filepath         - db dataset filepath will be converted to the same place as the xyz but replaced with .db format
            db_dataset                  - initialize AtomsData to create the db dataset with the filepath specified and available properties
                                          this class now can add molecules all at once using ase loaded molecules environment
                                          and the properties_list for all molecules (the dictionary defined above)
                                          each molecule should have a corresponding properties dictionary
                                          with all the available properties in that dinctionary
        
        Returns:
            a db file path next to the xyz file it was converted from 
            db files can be accessed with SchNet using AtomsData
    '''

    molecules = read(XYZDATASET_FILEPATH, index=':')
    
    number_properties = len(AVAILABLE_PROPERTIES)

    properties_list = []

    for mol in molecules:

        properties = [np.array([float(list(mol.info.keys())[each_prop])], dtype=np.float32) for each_prop in range(number_properties)]

        properties_list.append({AVAILABLE_PROPERTIES[each_prop]: properties[each_prop] for each_prop in range(number_properties)})
    
    db_dataset_filepath = XYZDATASET_FILEPATH.replace('.xyz','.db')
    db_dataset = AtomsData(db_dataset_filepath, available_properties=AVAILABLE_PROPERTIES)
    db_dataset.add_systems(molecules,properties_list)


def convert_xyz2mol(XYZ_FILENAME,MOL_FILENAME,BACKWARDS=False):
    '''
    uses os access to obabel to convert xyz files to mol files
    can also do it backwards

        Args:
            XYZ_FILENAME        - the filepath of the xyz molecule 
            MOL_FILENAME        - the filepath of where the mol should be converted to
            BACKWARDS           - boolean to convert from mol to xyz (False), or from xyz to mol (True)
        Returns
            converted mol file from xyz in the same directory (with same name but .mol)
    '''

    if BACKWARDS == False:
        os.system('obabel ' + XYZ_FILENAME +' -O ' + MOL_FILENAME)
    else:
        os.system('obabel ' + MOL_FILENAME +' -O ' +  XYZ_FILENAME)


def neighbortest_linearembanalogy_meandiff(REACTANTEMBS_FILEPATH,PRODUCTEMBS_FILEPATH,N_FEATURES,STACK_QM9,QM9EMBS_FILEPATH,LA=[]):
    '''
    Calculates the linear embedding analogy vector between reactant and product states as the mean
    of the difference between the embeddings of these two states.
    We assume that adding this vector produces the matching product embedding

        Args:
            REACTANT_EMBSFILEPATH        - the embeddings of the reactant target atoms
            PRODUCT_EMBSFILEPATH         - the embeddings of the product target atoms
        
        Process:
            reactants_embs               - reading in the reactants embs filepath as an array (matrix)
            product_embs                 - reading in the products embs filepath as an array (matrix) 
            diff_matrix_prodreact        - the difference matrix between products and reactants
            linear_analogy_vector        - the mean row of the difference matrix, which will act as the 
                                           linear analogy vector for the entire set (and possible other sets as well)
            transformed_embs             - reactant embs + linear_analogy_vector
            match_neighbortest           - method that checks if transformed embeddings is next to true product embedding
                                           meaning the linear analogy made sense in the embedding space. It returns the
                                           percentage of matches
            percent_match                - the percentage of matched neighbors embeddings, where the transformed embedding was nearest 
                                           to the true product embedding, in other words, percentage of reactions where linear anaolgy worked
            neighbor_idx_list            - the neighbor index of the matched embedding, 
                                           since the first transformed embedding should match the first product embedding
                                           and so does the second and third, fourth, and so on.... 
                                           then this list will inform you which neighbors matched and which did not because
                                           it should go in order 0,1,2,3,4,5,...N_molecules-1
                                           so if the order messes up that is a hint which reactions the linear analogy did NOT work
        Returns:
            linear_analogy_vector        - see above
            transformed_embs             - see above
            percent_match                - see above
            neighbor_idx_list            - see above

    '''
    
    reactant_embs = np.genfromtxt(REACTANTEMBS_FILEPATH,delimiter=',',encoding='utf-8-sig',skip_header=0)
    product_embs = np.genfromtxt(PRODUCTEMBS_FILEPATH,delimiter=',',encoding='utf-8-sig',skip_header=0)

    if LA == []:
        diff_matrix_prodreact = product_embs[:,:N_FEATURES] - reactant_embs[:,:N_FEATURES]
    
        linear_analogy_vector = np.mean(diff_matrix_prodreact, axis = 0)
    else:
        linear_analogy_vector = LA

    transformed_embs = reactant_embs[:,:N_FEATURES] + linear_analogy_vector

    percent_match, neighbors_idx_list = matching_neighbortest(transformed_embs,product_embs,STACK_QM9,QM9EMBS_FILEPATH,N_FEATURES)

    
    reactprod_embsstacked  = np.vstack((reactant_embs[:,:N_FEATURES],product_embs[:,:N_FEATURES]))
    return linear_analogy_vector, transformed_embs, percent_match, neighbors_idx_list, reactprod_embsstacked


def matching_neighbortest(TRANSFORMED_EMBS,PRODUCT_EMBS,STACK_QM9,QM9EMBS_FILEPATH,N_FEATURES):
    '''
    checks if the neighbor of the transformed embedding (by linear analogy as described above, using the mean difference vector)
    is the true embedding of the product. If it is then it is a match in the percentages of matches found that are neighbors with their
    true product embedding

        Args:
            TRANSFORMED_EMBS        - embeddings of the transformed product vis linead analogy
                                      as described above
            PRODUCT_EMBS            - embeddings of the true products, we will stack QM9 embs on top of this
                                      to make the task even more difficult, i.e the nearest even including vectors from QM9 embs
                                      of the same layer, would there still be a match? i.e putting the product_embs with the qm9_embs
                                      to see if confusion happens in the linear analogy!
            STACK_QM9               - if True, stacks qm9 embeddings (should be same layer) on top of the product embeddings to ensure 
                                      that the neighbor test involves those embeddings as well 
            QM9EMBS_FILEPATH        - filepath that contains QM9 test embeddings of the same layer as the ones that the product embeddings
                                      were extracted for
            N_FEATURES              - number of features in each embedding
            N_MOLECULES             - number of molecules to test

        Process:
            n_trans_embs            - number of transformed embeddings to perform neighbor test on
            qm9_embs                - if STACK_QM9 = True, then we load the qm9 embeddings
                                      as a numpy array from the QM9EBS_FILEPATH
            number_embs             - total number of embs in the stacked/non-stacked product embs 
                                      which will be used to check which is the nearest neighbor to the transformed embs
            percent_match           - a scalar that will compute the percent transformed that matched with their
                                      true product embedding, i.e had the nearest neighbor closest to the true embedding 
            neighbor_idx_list       - a list that will hold the neighbor indices, if the neighbor idx matches the product embedding idx
                                      then there is a match, a 100% match looks like [0,1,2,3,4,...] because all neighbor idx are matching 
                                      their product embedding, a mismatch in this list will point to where the linear analogy failed
            each_trans_emb          - an integer running through all the transformed embeddings
            distances_to_allembs    - list containing all the distances to other embs in product embs
                                      for each transformed embs
            neighbor_emb            - the index of the minumum distance is the neighbor index
                                      if this equals the true product embeddings, 
                                      which are stacked below qm9 embeddings, then there is a match 
                                      and the linear analogy worked


        Returns:
            percent_match           - see above
            neighbor_idx_lise       - see above
    '''

    n_trans_embs = len(PRODUCT_EMBS)
    
    if STACK_QM9 == True:
        qm9_embs = np.genfromtxt(QM9EMBS_FILEPATH,delimiter=',',skip_header=1)
        PRODUCT_EMBS = np.vstack((PRODUCT_EMBS[:,:N_FEATURES],qm9_embs[:,:N_FEATURES]))
    
    number_embs = len(PRODUCT_EMBS)

    percent_match = 0
    neighbor_idx_list = []
    for each_trans_emb in range(n_trans_embs):

        distances_to_allembs = [np.linalg.norm(TRANSFORMED_EMBS[each_trans_emb,0:N_FEATURES] - PRODUCT_EMBS[each_emb,0:N_FEATURES]) for each_emb in range(number_embs)]
        
        neighbor_emb = np.argmin(distances_to_allembs)
        
        neighbor_idx_list.append(neighbor_emb)

        if neighbor_emb == each_trans_emb:
            percent_match = percent_match + 1
    
    percent_match = (percent_match/(n_trans_embs))*100
    
    return percent_match, neighbor_idx_list


def connect_reactprod_scatters2vectors(DATA_FILEPATH,SAVE_FILEPATH):
    '''
    this transforms our scatter plot of embeddings (preferably in PC basis)
    to reactants/products vectors that can be plotted with gnuplot (x_r,y_r,x_p,y_p)
    x_r, y_r are found in the first half of the scatter data
    x_p, y_p are found int he second half of the scatter data

        Args:
            DATA_FILEPATH       - the filepath where the scatter pca of reactant + product embeddings
                                  can be found
            SAVE_FILEPATH       - where to save the resulting vector-to-vector data
                                  (x_r,y_r,x_p,y_p)

        Process:
            data                - np.array of the data in the DATA_FILEPATH (np.genfromtxt)
            vector_data         - the np.array that will hold the vector data described above
            each_data           - an integer running over only half of the data, 
                                  it runs over the first half and second half simultaneously
                                  so that it is able fill the vector data made up up
                                  of reactant and product embeddings scatter
                                  
        Returns:
            saved vector-to-vector data in the specified csv filepath
    '''
    #generate data
    data = np.genfromtxt(DATA_FILEPATH,delimiter=',',encoding='utf-8-sig',skip_header=1)

    vector_data = np.zeros((int(len(data)/2),4))
    #Assuming hald the data is init, half is final
    for each_data in range(0,int(len(data)/2)):
        vector_data[each_data,0:2] = data[each_data,0:2]
        vector_data[each_data,2:4] = data[each_data+int(len(data)/2),0:2] - vector_data[each_data,0:2]

    np.savetxt(SAVE_FILEPATH,vector_data,delimiter=',')

def connect_reactprod_scatters2vectors_fullvector(DATA_FILEPATH,SAVE_FILEPATH):
    '''
    this transforms our scatter plot of embeddings (preferably in PC basis)
    to reactants/products vectors that can be plotted with gnuplot (x_r,y_r,x_p,y_p)
    x_r, y_r are found in the first half of the scatter data
    x_p, y_p are found int he second half of the scatter data

        Args:
            DATA_FILEPATH       - the filepath where the scatter pca of reactant + product embeddings
                                  can be found
            SAVE_FILEPATH       - where to save the resulting vector-to-vector data
                                  (x_r,y_r,x_p,y_p)

        Process:
            data                - np.array of the data in the DATA_FILEPATH (np.genfromtxt)
            vector_data         - the np.array that will hold the vector data described above
            each_data           - an integer running over only half of the data, 
                                  it runs over the first half and second half simultaneously
                                  so that it is able fill the vector data made up up
                                  of reactant and product embeddings scatter
                                  
        Returns:
            saved vector-to-vector data in the specified csv filepath
    '''
    #generate data
    data = np.genfromtxt(DATA_FILEPATH,delimiter=',',encoding='utf-8-sig',skip_header=1)

    vector_data = np.zeros((int(len(data)/2),256))
    #Assuming hald the data is init, half is final
    for each_data in range(0,int(len(data)/2)):
        vector_data[each_data,0:128] = data[each_data,0:128]
        vector_data[each_data,128:256] = data[each_data+int(len(data)/2),0:128] - vector_data[each_data,0:128]

    np.savetxt(SAVE_FILEPATH,vector_data,delimiter=',')


def find_if_reactioncenter_exists_exactconfig(CENTER_STRUCTURES,CENTER_LABELS,MOLECULE_LABELS):
    '''
    A central tool, used multiple times, to find out if atom center structures sepcified in the input
    for reaction centers, leaving centers, and centers of adding/removal of double bonds
    are found with the specific number they are allowed to be found in.
    The code uses all the possible labels that represent the center structure found by searching for the center structure's
    Possible labels in depth 5 or more label file

    Args:
        CENTER_STRUCTURES           - follows the format of the input for specifying structures and number of of that structure we are looking for:
                                        [['center_structure1',#1],['center_stucture2',#2]]
        found_toomany_rxnycenter    - a bool flag that will be used to break out of the loop of finding center labels matches in the molecule,
                                      if any center has too many (more than specified for it) matched atom-labels, the flag allows to break to the next molecule
                                      otherwise we would just be breaking out of the loop finding the label matches amongst all the possible labels for the center structure 
        does_eachcenter_exist       - a list of booleans ('yes' or 'no') that will hold the bool value for IF the center structure has been found to exist (via center labels) 
                                      in the exact number specified for it in the input
        each_center                 - running through each center as an integer because it cannot be easily accessed otherwise
        matchedlabels_center        - list that holds the matched labels for the atom centers that have been found to match any of the possible
                                      centers labels derived from the center structure given in the input
                                      this list will be checked each time a new match is to be added, if the number of unique matches, surpasses the number
                                      specified for the input center structure the algorithm breaks for that whole molecule and goes to the next
        center_indices              - list that will hold all the indices of the atom centers that have the matched labels
        each_center                 - integer running through the center structures (otherwise it is not as easy to access directly)
        label                       - integer running through the labels possible for each center structure (also easiest to access this way)
        MOLECULE_LABELS             - the N-DEPTH labels of each atom in the molecule, used to find matches to the center structure labels
    '''

    #A flag that will be used to break out of the next loop, if one reaction 
    #center has too many matched atom-labels then break to the next molecule
    found_toomany_rxncenter = False
                
    #a list that will hold booleans, 'yes' or 'no' depending
    #on IF the label for the atom reaction centers
    #are found in the labels of this molecule
    does_eachcenter_exist = []

    #IF the reaction center's labels have been found in the molecule
    #the labels for it will be held in this list
    #the number of labels collected in this list will be checked every time we add a label to it
    #the number must never surpass that which is allowed for the reaction structure specified next 
    #to each reaction center structure in the input
    #NOTE the list may be a variety of labels representing a variety of centers that have the same allowed structure
    matchedlabels_centers = []

    #All idxs in the molecule where ALL the structure centers have been found will be collected here in ONE list
    all_centers_idxs = []


    #for each reaction center, you need to find a match for it in the envlabels, 
    #the number of matches in the molecule must also be equal to the number specified next to the structure
    for center_i in range(len(CENTER_STRUCTURES)):
    

        #list that holds the matched labels for this center
        #this will be checked each time a new atom-label is added as a possible center
        #to ensure that the number of matched labels for this center does not exceed the number specified in the input
        matchedlabel_thiscenter = []

        #RUN through all the possible labels for the specified reaction center
        for label_i in range(len(CENTER_LABELS[center_i])): 
            
            #IF one of the labels for the structure is found in envlabels of molecule
            if CENTER_LABELS[center_i][label_i] in MOLECULE_LABELS.tolist():
                
                #IF we don't have more than the allowed number of matches specified for that 
                #reaction center than we can allow the match, otherwise we have to stop and go to the
                #next  molecule, as this molecule contains too many available reaction centers
                if len(matchedlabel_thiscenter) > CENTER_STRUCTURES[center_i][1]:  #checking that there are not too many reaction centers

                    found_toomany_rxncenter = True
                    #reset match label to be empty
                    matchedlabels_centers = []
                    matchedlabel_thiscenter = []
                    #start for the next atom, this atom has too many reaction centers
                    break
                else:
                    matchedlabel_thiscenter.append(CENTER_LABELS[center_i][label_i])
                    #one of the possible labels matched for this reaction center
                    matchedlabels_centers.append(CENTER_LABELS[center_i][label_i])

        #break out of the molecule completely, if you found any ONY center that has too many labels
        if found_toomany_rxncenter:
            break

        #checking that  the reaction center has exactly the number of matches acceptable
        #NOT less for example (as that can happen, the above conditions only eliminate above number)
        #NOTE below even works for the symmetric case, where same label is found on two atoms, long as the number of 
        #counts of atoms labels symmetric exactly/non-symmetric equals the allowed for that center then it is a "yes" found center
        if np.count_nonzero(MOLECULE_LABELS[np.isin(MOLECULE_LABELS, matchedlabel_thiscenter)]) == CENTER_STRUCTURES[center_i][1]:
            #Bool: This reaction center exists! WITH exactly the number of allowed matches (of any of possible labels for this structure)
            #will be used later to ensure that all our reaction centers exist as specified
            does_eachcenter_exist.append('yes')

            #collect indices where the matched labels occur for each center center
            center_indices = np.where(np.isin(MOLECULE_LABELS, matchedlabel_thiscenter))[0].tolist()


            for center_idx in center_indices:
                all_centers_idxs.append(center_idx)

        else:
            #Bool: This reaction center does not exist or meet the conditions
            #used to find out which reaction centers don't exist 
            #used to find out which reaction centers don't exist 
            #or don't have the number of matched labels
            #This means that the whole reaction center is NOT RIGHT
            does_eachcenter_exist.append('no')  

        
    return does_eachcenter_exist, matchedlabels_centers, all_centers_idxs

def find_leavingatoms(CENTER_STRUCTURES,CENTER_LABELS,MOLECULE_LABELS):
    '''
    This code helps to find leaving atoms around a specified reaction center that has already been found uniquely
    ONLY here finding uniqueness is not a necessity, this code finds a number of specified atoms around the reaction cetner 

    Args:
        CENTER_STRUCTURES           - the atom-environment structures to be looking for, 
                                      they also specify how many to find but they 
                                      don't have to be found uniquly on the molecule 
                                      (as with the reaction center)
        CENTER_LABELS               - all the possible N-depth atom-environment labels for each
                                      specified center structure 
        MOLECULE_LABELS             - labels of atom environments of the molecule 
                                      (at a certain depth, obtained from the autolabel code)
        
        
                                      
    Process:
        does_center_exist           - finds if each atom env structure is found on the molecule (ONLY ONCE! see note below)
                                      IF THIS IS DEFINED PROPERLY BASED AROUND THE UNIQUE REACTION CENTER
                                      THEN YOU SHOULD RECEIVE NO ERRORS
                                      must be done right in order to remove the correct leaving atoms around the reaction center
        indices_allcenters_found    - indices of the molecule at which the center was found, used for locating purposes, 
                                      helps to locate the atoms to remove or to change bond order of in later parts of the code
        center_i                    - running through the specified atom centers to remove
        matched_label               - will hold the label found for the structure to be removed 
                                      NOTE there should ONLY be one matched label per structure specified, that is because there is ONLY 
                                      ONE REACTION CENTER ALLOWED, and therefore if specified right, there will only be one of this structure
                                      which is based around the reaction center (an atom to leave)
        label                       - the possible integer labels of each structure
                                     (found previously from searching the structure in the depth N atom env label file )
        indices_matched_label       - indices of the molecule where the matched label has been found
        numberidxs_tohold           - some instances, due to symmetry, we find a structure, 
                                      but only want to remove one atom for ex. the H in H2C, 
    '''
    does_center_exist = []
    indices_allcenters_found = []
    labels_allcenters_found = []

    for center_i in range(len(CENTER_STRUCTURES)):
        matched_label = None

        #First check that any one and ONLY one of the possible labels for a leaving atom is found
        #BUT you do not necessarily need to check that it is found only once in the molecule
        for label in CENTER_LABELS[center_i]:  
            if label in MOLECULE_LABELS.tolist():
                if matched_label is not None:  # Another match found, so not unique
                    matched_label = None
                    break
            
                else:
                    matched_label = label
                    labels_allcenters_found.append(matched_label)

        if matched_label is not None:
            #For the unique label found out of all the possibilities of labels for the substructure at this depth
            #find the indice in the molecule that match this label
            indices_matched_label = np.where(MOLECULE_LABELS == matched_label)[0]
            numberidxs_tohold = CENTER_STRUCTURES[center_i][1]

            #FOR SPECIAL CASES WHEN H's beside each other have to be avoided
#            indices_matched_label = indices_matched_label[::2]

            for each_idx in range(numberidxs_tohold):
                indices_allcenters_found.append(indices_matched_label[each_idx])

            does_center_exist.append('yes')
        else:
            does_center_exist.append('no')

    return does_center_exist, labels_allcenters_found, indices_allcenters_found


def find_bondorderchange_atoms(CENTER_STRUCTURES,CENTER_LABELS,MOLECULE_LABELS):
    '''
    Finds atom indices to change bond order of as part of the reaction process

    NOTE WHILE THIS FUNCTION LOOKS THE SAME AS THE ABOVE TWO, IT IS DIFFERENT FROM BOTH IN VERY SUBTLE WAYS
    IMPORTANT TO KEEP AS DISTINCT FUNCTION AND NOT COMBINE!!! HERE, ONE CENTER STRUCTURE
    IS ALLOWED TO HAVE MULTIPLE LABELS, BECAUSE BONDS ARE FOUND ACROSS ATOMS AND THE CENTER STRUCTURE
    MAY BE SYMMETRIC FOR THAT PURPOSES. IF YOU ONLY ALLOW ONE MATCHED LABEL (AS FOR LG) IT WILL NOT WORK
    HOWEVER, SIMILAR TO LG CODE BECAUSE UNIQUENESS IS NOT NECESSARY THIS IS BECAUSE YOU MIGHT FIND TOO MANY OF ONE CENTER BECAUSE OF SYMMETRY (THIS HAPPENS IN CYCLOPROPANOL FOR EXAMPLE)
    AND STILL BE ABLE TO FORM A BOND ACROSS.... JUST DONT FUCKING DELETE IT

    Args:
        CENTER_STRUCTURES           - the atom-environment structures to be looking for, 
                                      they also specify how many to find but they 
                                      don't have to be found uniquly on the molecule 
                                      (as with the reaction center)
        CENTER_LABELS               - all the possible N-depth atom-environment labels for each
                                      specified center structure 
        MOLECULE_LABELS             - labels of atom environments of the molecule 
                                      (at a certain depth, obtained from the autolabel code)
        
    Process:
        indices_changebond          - holds the indices of the molecule at which the bond order has to change (found through searching for structure label)
        allmatched_labels           - holds the matched labels of the center structure to add/rem bond 
        center_i                    - running through the specified atom centers to remove
        matched_label               - will hold the label found for the structure to be changed bond order of
                                      NOTE there should ONLY be one matched label per structure specified, that is because there is ONLY 
                                      ONE REACTION CENTER ALLOWED, and therefore if specified right, there will only be one of this structure
                                      which is based around the reaction center (an atom to leave)
        label                       - the possible integer labels of each structure
                                     (found previously from searching the structure in the depth N atom env label file )
        indices_matchedlabel        - indices of molecule at which the center's label have been matched
        indices_changebond          - indices of the molecule where the matched label has been found
    
    '''

    indices_changebond = []
    allmatched_labels = []

    for center_i in range(len(CENTER_STRUCTURES)):
        matched_label = None

        #First check that any one and ONLY one of the possible labels for a leaving atom is found
        #BUT you do not necessarily need to check that it is found only once in the molecule
        for label in CENTER_LABELS[center_i]:  
            if label in MOLECULE_LABELS.tolist():
                if len(allmatched_labels) > CENTER_STRUCTURES[center_i][1]:  # ENOUGH matched found, stop
                    matched_label = None
                    break
                else:
                    matched_label = label
                    allmatched_labels.append(matched_label)

        if matched_label is not None:
            #For the unique label found out of all the possibilities of labels for the substructure at this depth
            #find the indice in the molecule that match this label
            indices_matched_label = np.where(MOLECULE_LABELS == matched_label)[0]

            if len(indices_matched_label) > 1:
                indices_changebond.append(indices_matched_label)
            elif len(indices_matched_label) == 1: 
                indices_changebond.append(indices_matched_label[0])


    return indices_changebond