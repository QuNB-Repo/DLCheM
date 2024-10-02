
from schnetpack.datasets import QM9
from schnetpack import AtomsData
import pandas as pd
import csv
import numpy as np

from runrxns.utils import utils, utilssym

from rdkit import Chem
from rdkit.Chem import AllChem

import warnings

warnings.filterwarnings("ignore")

'''
Last Updated: 2024-03-19

'''

class RemoveBuild_LGFG_SingleBuildSite():
    '''
    transforms reactants to products by cutting atoms at the reaction center (leaving group), 
    then attaching incoming functional groups, and maybe adds/removes a double bond at the reaction center
    the ultimate result is to change a reactant form of a molecule to its product form mimicking a reaction
    
        functions in the class

            __init__             - initializes the control and build variables
            Write_InitXYZString  - writes the reactant molecule from db file to a string, 
                                   if appropriate reaction center is found on reactant molecule
            RemoveAtoms_LG       - function that removes atoms at the reaction center from an xyz file
                                   returns the resulting xyz without any optimization for analysis of embeddings
                                   of this intermediatte step (before adding the functional group),
                                   a.k.a the leaving group step
            Build_FG             - function that attaches the functional group to the reaction center 
                                   and optimizes it using RDKit's MMFF94, 
                                   NOTE requires conversion to mol file with obabel because functional groups
                                   have to be tacked on to the correct place (build_index) 
            RemoveBuild_LGFG     - the function that works other functions above to find and remove atom labels 
                                   associated with an reaction center and then build a functional group
                                   instead of the removed leaving group
    '''

    def __init__(self,DATASET_FILEPATH,LABELS_FILEPATH,DEPTH,RXN_CENTERS,LEAVING_ATOMS,ATTACH_GROUP,ATTACH_GROUP_BONDS,ATTACH_GROUP_NUMB_ATOMSBONDS,ADD_BOND_TO_RXNSITE,REMOVE_BOND_FROM_RXNSITE,MOLECULES_RANGE,AVAILABLE_PROPERTIES,QM9_BOOL,REACT_FILEPATH,REMOVELG_FILEPATH,PROD_FILEPATH,SKIP,SKIP_SYMMETRIC_ATOMS):
        '''
        Initializes the control  & build variables:

            control:

                DATASET_FILEPEATH           - the db file that has been labelled with autolabel and which autolabel file is located somewhere known
                LABELS_FILEPATH             - the filepath of the atomic environment labels at sufficient depth (usually a depth of 5 is enough to pick up most reaction centers)
                RXN_CENTERS                 - the structure labels that identify the atomic centers of reaction, where the leaving group leaves, and where
                                              the incoming FG attaches. It is possible to have multiple atoms that are part of the reaction site (ex. diels alder)
                                              the relevant parts of the atomic center structure have to be identified by the label and the rest left open for freedom and use regular expressions where possible
                                              SEE fgtransform file that calls this code to read about how to define a reaction center structure (or pick one out from a preset dictionary using a reaction center name) 
                LEAVING_ATOMS               - specifies which atoms (using local structures labels at certain depth) need to be removed (from or not-from) the reaction center structure
                                              the algorithm will use the structure label to find the exact
                                              atoms to remove based on the specified structure around them  
                                              SEE fgtransform (which calls this code) to read about how to define leaving group atom structure labels
                                              or pick one out from a preset dictionary using leaving group name                

                ATTACH_GROUP                - this is mol string information about the functional group to build
                                              it is the line of fg to build, usually all atoms at coordinates 0,0,0 but will be optimized
                                              the bonding is what helps the force field and that is defined separately below. And example of the 
                                              this part would be 

                                              for and OH functional group:
                                              '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n'

                                              it is possible to specify a none FG, if there is no FG to build:
                                              '' 

                                              this tells the algorithm to add nothing but optimize the 
                                              leaving-group-removed structure (important for oxidation of alcohols, as sometimes 
                                              there is no functional group to build, there is just functional group leaving)
                                              NOTE that even without a building functional, you should still specify a build label
                                              this helps the embedding extraction code pick out which atom's embedding to extract
                                              to analysze
                                              NOTE in the next update make a separate 
                                              list object to tag all the atom-indices 
                                              you want to analyze with embeddings
                ATTACH_GROUP_BONDS          - a template for the bonding section of the mol file that will replace alphabetical symbols
                                              with the correct connections
                                            
                                              ex.
                                              fg_OH_bond = ' a f  1  0\n b a  1  0\n'

                                              the code knows that 'a' means index 1, 'b' means 2, and so on
                                              with the exception of 'f' which is replaced with the bond_index of the functional group
                                              found by using the atombuild_label
                                              this is kept as a template copy and the Build_FG method copies it 
                                              replace the alphabets with the correct connecting indices
                ATTACH_GROUP_NUMB_
                ATOMSBONDS                  - this specifies to the algorithm how many atoms and bonds this functional group adds 
                                              so it can adjust the mol file accordingly

                                              ex.
                                              [2,2]    - add two atoms, add two bonds
                ADD_BOND_TO_RXNSITE         - specifies which atomic centers (using structure labels) to add a double bond around two centers if required
                                              some reactions require this (ex. oxidation)
                REMOVE_BOND_FROM_RXNSITE    - specifies which atomic center (using structure labels) to remove a double bond around two centers if required
                                              some reactions requires this (ex. diels alder)
                MOLECULES_RANGE             - range of molecules from the db file to scan for atomic environments and thus allow 
                                              through the algorithm for transformation (removing atomic environments, adding functional groups)
                                              
                                              NOTE this range needs to be labeled by autolabel in the label_filepath!
                
                AVAILABLE_PROPERTIES        - this is required for making the db files from xyz files,
                                              db files comes with associated properties, 
                                              for now this is a filler! and can be anything! just keep note of what name it is so the db can be accessed later
                QM9_BOOL                    - a boolean depending on whether qm9.db is being used because to load all the properties of that one
                                              we use the QM9 method otherwise the AtomsData method from schnetpack
                REACT_FILEPATH              - the filepath to create the data bse of reactants file
                REMOVELG_FILEPATH           - the filepath to create the database of reactants w/ leaving group removed WITHOUT optimization
                PROD_FILEPATH               - the filepath to create the database of corresponding products with optimization
                label                       - the atom labels that define atomic structure at a particular depth. Will be used by the algorithm
                                              to find appropriate reaction sites, leaving groups, and attach functional group to the right place      
                atomlabels_rxncenters       - list of lists that will hold all depth N structure integer labels that contain the substructures specified by each reaction site   
                                              each reaction site specified gets a list of integers of depth N labels that are appropriate for it                                         
                atomlabels_leavinggroup     - list of lists that will hold all depth N structure integer labels that contain the substructures specified by each leaving atom   
                                              each leaving atom specified gets a list of integers of depth N labels that are appropriate for it 
                atomlabels_addbond          - list of lists that will hold all depth N structure integer labels that contain the substructures specified by each atomic site to add bond to  
                                              each atom around the bond site specified gets a list of integers of depth N labels that are appropriate for it 
                atomlabels_removebond       - list of lists that will hold all depth N structure integer labels that contain the substructures specified by each atomic site to remove a bond from  
                                              each atom around the bond site specified gets a list of integers of depth N labels that are appropriate for it 
                each_atom                   - an integer running through either the number of reaction sites/leaving atoms/bond centers
                unique_labels               - filtered integer labels from the label depth file that represent structure that contain the substructure of reaction sites/leaving atoms/bond centers

            build:
                react_dataset_xyzfile       - open the init file in write mode
                removelg_dataset_xyzfile    - open the rem file in write mode
                prod_dataset_xyzfile        - open the trans file in write mode
                all_build_indices           - this is important to tag the build_index in a list so that
                                              the extractembeddings file can extract embeddings of
                                              the build atom 
                                              (where the leaving group is detaching and attacking group attacking, for example)
                all_build_initindices       - after removal of the leaving, the tagged build index may be different
                                              than before the leaving group, this ensures a tag before the leaving group
                                              for the embeddings of the reactant atom to be extracted
                all_targetextraction_
                initindices                 - NOTE  NOT YET IMPLEMENTED, in case you want to examine embeddings of other sites than the reaction sites
                                              this list of lists will hold all indices of the target site for each reactant molecule
                self.all_targetextraction_
                indices                     - NOTE NOT YET IMPLEMENTED, in case you want to examine embeddings of other sites than the reaction sites,
                                              this list of lists will hold all indices of the target site for each prod/removelg molecule                             
        '''

        #control variables
        self.DATASET_FILEPEATH = DATASET_FILEPATH
        self.LABELS_FILEPATH = LABELS_FILEPATH
        self.DEPTH = DEPTH
        self.RXN_CENTERS  = RXN_CENTERS
        self.LEAVING_ATOMS = LEAVING_ATOMS
        self.ATTACH_GROUP = ATTACH_GROUP
        self.ATTACH_GROUP_BONDS = ATTACH_GROUP_BONDS
        self.ATTACH_GROUP_NUMB_ATOMSBONDS = ATTACH_GROUP_NUMB_ATOMSBONDS
        self.REMOVE_BOND_FROM_RXNSITE = REMOVE_BOND_FROM_RXNSITE
        self.ADD_BOND_TO_RXNSITE = ADD_BOND_TO_RXNSITE
        self.MOLECULES_RANGE = MOLECULES_RANGE
        self.AVAILABLE_PROPERTIES = AVAILABLE_PROPERTIES
        self.QM9_BOOL = QM9_BOOL
        self.REACT_FILEPATH = REACT_FILEPATH
        self.REMOVELG_FILEPATH = REMOVELG_FILEPATH
        self.PROD_FILEPATH = PROD_FILEPATH
        self.SKIP = SKIP
        self.SKIP_SYMMETRIC_ATOMS = SKIP_SYMMETRIC_ATOMS

        #load the labels at a sufficient depth
        self.labels = pd.read_csv(self.LABELS_FILEPATH)

        self.atomlabels_rxncenters = []
        for each_atom in range(len(self.RXN_CENTERS)):
            unique_labels = np.unique(self.labels[self.labels['key'+str(self.DEPTH)].str.match(self.RXN_CENTERS[each_atom][0])]['labeldepth'+str(self.DEPTH)].values).tolist()
            self.atomlabels_rxncenters.append(unique_labels)

        self.atomlabels_leavinggroup = []
        if self.LEAVING_ATOMS is not None:
            for each_atom in range(len(self.LEAVING_ATOMS)):
                unique_labels = np.unique(self.labels[self.labels['key'+str(self.DEPTH)].str.match(self.LEAVING_ATOMS[each_atom][0])]['labeldepth'+str(self.DEPTH)].values).tolist()
                self.atomlabels_leavinggroup.append(unique_labels)

        self.atomlabels_removebond = []
        if self.REMOVE_BOND_FROM_RXNSITE is not None:
            for each_atom in range(len(self.REMOVE_BOND_FROM_RXNSITE)):
                unique_labels = np.unique(self.labels[self.labels['key'+str(self.DEPTH)].str.match(self.REMOVE_BOND_FROM_RXNSITE[each_atom][0])]['labeldepth'+str(self.DEPTH)].values).tolist()
                self.atomlabels_removebond.append(unique_labels)

        self.atomlabels_addbond = []
        if self.ADD_BOND_TO_RXNSITE is not None:
            for each_atom in range(len(self.ADD_BOND_TO_RXNSITE)):
                unique_labels = np.unique(self.labels[self.labels['key'+str(self.DEPTH)].str.match(self.ADD_BOND_TO_RXNSITE[each_atom][0])]['labeldepth'+str(self.DEPTH)].values).tolist()
                self.atomlabels_addbond.append(unique_labels)

        if len(self.atomlabels_rxncenters) != len(self.RXN_CENTERS):
            raise KeyError("Could not find one or more of the reaction center structures specified, please check how you are specifying reaction center structure")
        if self.LEAVING_ATOMS is not None:
            if len(self.atomlabels_leavinggroup) != len(self.LEAVING_ATOMS):
                raise KeyError("Could not find one or more of the leaving atoms structures specified, please check how you are specifying your leaving atom structures")
        if self.ADD_BOND_TO_RXNSITE is not None:
            if len(self.atomlabels_addbond) != len(self.ADD_BOND_TO_RXNSITE):
                raise KeyError("Could not find one or more of the add-bond atoms structures specified, please check how you are specifying your add-bond atom structures")
        if self.REMOVE_BOND_FROM_RXNSITE is not None:
            if len(self.atomlabels_removebond) != len(self.REMOVE_BOND_FROM_RXNSITE):
                raise KeyError("Could not find one or more of the remove-bond atoms structures specified, please check how you are specifying your remove-bond atom structures")

        #build variables
        #these filepaths will be the same place where xyz will be converted to the db file
        self.react_dataset_xyzfile = open(self.REACT_FILEPATH, mode='a')
        self.removelg_dataset_xyzfile = open(self.REMOVELG_FILEPATH,mode='a')
        self.prod_dataset_xyzfile = open(self.PROD_FILEPATH, mode='a')

        #save all_build_indices for each molecule,
        #there will be ONLY 1 per molecule, 
        #this is so that we can use it to find 
        #the embedding of the target build atom
        #for init_ it has to be found by label,
        # as build_label may be different before/after removal
        self.all_reactcenterindices = []
        self.all_prodcenterindices = []
        
        #You can set your own target, for embedding extraction
        #it doesn't have to be the build_index of each molecule, for later
        #NEXT UPDATE
        self.all_targetextraction_initindices = []
        self.all_targetextraction_indices = []



    def Write_InitXYZString(self,atomic_numbers,numb_atoms,xyz_positions):
        '''
        Writes the reactant if the atomic environment has been found in the db file

            react_xyz_string         - the xyz string that will be written based on xyz
                                      information about the molecules
            numb_atom               - number of atoms, info comes from the db file
            xyz_positions           - xyz coordinates, info comes from the db file
            each_atom               - an integer (0 -> N_atoms-1) for each atom in molecule
            atomic_numbers          - atomic numbers list, info comes from the db file
        ''' 
        
        react_xyz_string = ''
        react_xyz_string = react_xyz_string + str(numb_atoms) + '\n0.0000 \n'
        for each_atom in range(numb_atoms):
            if atomic_numbers[each_atom] == 1:
                react_xyz_string = react_xyz_string + 'H' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 6:
                react_xyz_string = react_xyz_string + 'C' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 7:
                react_xyz_string = react_xyz_string + 'N' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 8:
                react_xyz_string = react_xyz_string + 'O' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 9:
                react_xyz_string = react_xyz_string + 'F' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
        return react_xyz_string

    def RemoveAtoms_LG(self,indices_remove):
        '''
        function that removes the leaving group from the molecule 
        (as specified by atomlabels_remove)

            removelg_xyz_string           - the resulting xyz_string which will have the leaving group removed
            count_xyzline               - initial counter for how many xyz lines have been passed
            this_line_index             - indices of lines in the xyz file
            this_line                   - lines in the xyz file
            react_xyz_string             - the initial xyz string which needs to be modifed by this function 
                                          (by removing leaving group)
            indices_remove              - the indices of the atoms to remove 
                                          (obtained previously from using the atomlabels_remove found for this molecule)
            removelg_xyz_withemptylines   - the lines of the removed_xyz_string which may have
                                          some lines empty, now that the leaving group has been removed
            removelg_xyz_stringnonempty_
            lines_list                  - a list that has all the lines in a 
            removelg_xyz_string_noempty   - it is important for this to work 
                                          that the resulting string is not copied on top
                                          of removelg_xyz_string, we use this copy to join the lines_list 
                                          above with '\n' thus giving us the xyz with LG removed
                                          WITH NO EMPTY LINES

                                          IT WILL NOT WORK IF YOU USE ORIGINAL BECAUSE
                                          IT WILL THINK TO JOIN '\n' the lines from the the original as well
        '''

        removelg_xyz_string = ''
        count_xyzline = 0

        for this_line_index, this_line in enumerate(self.react_xyz_string.split('\n')):
            
            if this_line_index == 0:
                this_line = this_line.replace(this_line,'\n'+str(int(this_line)-len(indices_remove)) + '\n')
            
            if this_line_index > 1: 
                if count_xyzline in indices_remove:
                    this_line = this_line.replace(this_line,'')
                count_xyzline = count_xyzline + 1

            removelg_xyz_string = removelg_xyz_string + this_line + '\n'

        #remove empty lines
        removelg_xyz_withemptylines = removelg_xyz_string.split('\n')
        removelg_xyz_stringnonempty_lines = [line for line in removelg_xyz_withemptylines if line.strip()]
        # Join the non-empty lines back into a string
        removelg_xyz_string_noempty = '\n'.join(removelg_xyz_stringnonempty_lines)+'\n'

        return removelg_xyz_string_noempty

    def Build_FG(self,build_indices,number_atoms,indices_addbond,indices_removebond):
        '''
        This is the part that builds a fully specified functional group (FG_TO_BUILD,FG_BOND_TO_BUILD,ADD_ATOMS_BONDS)
        to the xyz file, it will have to do a temporary conversion to a mol file to be able to do the bond connections
        and then perform geometry optimization with MMFF94 using RDKit

            build_index                 - the index where to build the FG (obtained from ATOMLABELS_BUILD, ex. ['labeldepth2',32'])
            temp_xyz_filename           - a temporary filepath to write the xyz so it can be converted to mol
            temp_xyz_file               - opening the temporary file in write mode
            temp_mol_filename           - a temporary mol filename for os with obabel to convert from xyz to this mol filename
            utils.convert_xyz2mol       - uses obabel with os to convert temp xyz to mol file
            mol                         - the mol in RDKit, loaded from the mol_file
            mol_block                   - the mol_block (i.e string) using Chem.MolFromMolBlock
            FG_BOND_TO_BUILD            - this comes from input and specifies the bond template
                                          the bond template will be replaced with the actual bonding
                                          indices that will connect the FG with molecule. also needs to contain 
                                          bonding indices of the new FG
                                          There is a alphabet:index key that makes this work
                                            FG will always go at the end of the mol/xyz file:

                                                'a' - number_atoms + 1
                                                'b' - number_atoms + 2
                                                'c' - number_atoms + 3
                                                'd' - number_atoms + 4
                                                'z' - bond_index
                                         
                                           'a' tells where the first atom of the FG, 'b' the second atom of the FG and so on
                                           this comes in straight from the input 

                                           'z' is used to identify the bond_index where the FG will connect
                                            with the original molecule

                                          NOTE important NOT to update self.FG_BOND_TO_BUILD (input control), as that needs to
                                          stay a template for other rounds, instead save the replacements in "FG_BOND_TO_BUILD" NO update
            number_atoms                - number of atoms in the FG removed file, after the LG is removed
            each_line_index             - indexing each line in the mol file
            each_line                   - lines in the mol file
            number_bonds                - this is read off from the mol file at line 3, int(each_line[4:6]) is the bond number
                                          for the molecule with FG removed
            new_number_atoms            - comes from input control, first number in the ADD_ATOMS_BONDS list 
                                          [add_atoms,add_bonds]
            new_number_bonds            - the second number in the list above, tells the algorithm how many bonds to add
                                          mol file requires a listing of both number of atoms and number of bonds
            FG_TO_BUILD                 - comes from input, tells the algorithm which FG to build in mol format, see input control
            temp_out_molfilename        - a temporary scratch file for writing  optimized mol file so that it can be converted to xyz 
                                          with os and obabel
            temp_out_xyzfilename        - a temporary xyz file that holds the converted mol file 
            opt_xyz_file                - opening the temporary xyz file that has the optimized FINAL structure
            opt_prod_xyz_string        - writing the optimized final structure as a string and returning that
        '''

        build_indices  = (np.array(build_indices)+1).tolist() #to match the indexing of mol file bonding section which starts at 1

        temp_xyz_filename = 'temp.xyz'
        temp_xyz_file = open(temp_xyz_filename, mode='w')
        temp_xyz_file.write(self.removelg_xyz_string)
        temp_xyz_file.close()

        temp_mol_filename = 'temp.mol'
        utils.convert_xyz2mol(temp_xyz_filename,temp_mol_filename)

        mol = Chem.MolFromMolFile(temp_mol_filename,sanitize=False,removeHs=False)
#        mol = Chem.AddHs(mol)

        mol_block = Chem.MolToMolBlock(mol)

        fg_bond_to_build = self.ATTACH_GROUP_BONDS 

        if number_atoms + 1 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('a',' '+str(number_atoms+1))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('a',''+str(number_atoms+1))
        
        if number_atoms + 2 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('b',' '+str(number_atoms+2))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('b',''+str(number_atoms+2))
        
        if number_atoms + 3 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('c',' '+str(number_atoms+3))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('c',''+str(number_atoms+3))
        
        if number_atoms + 4 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('d',' '+str(number_atoms+4))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('d',''+str(number_atoms+4))
        if number_atoms + 5 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('e',' '+str(number_atoms+5))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('e',''+str(number_atoms+5))
        
        if number_atoms + 6 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('f',' '+str(number_atoms+6))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('f',''+str(number_atoms+6))
        
        if number_atoms + 7 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('g',' '+str(number_atoms+7))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('g',''+str(number_atoms+7))
        
        if number_atoms + 8 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('h',' '+str(number_atoms+8))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('h',''+str(number_atoms+8))
        
        if number_atoms + 9 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('i',' '+str(number_atoms+9))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('i',''+str(number_atoms+9))
        
        if number_atoms + 10 < 10:
            fg_bond_to_build = fg_bond_to_build.replace('j',' '+str(number_atoms+10))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('j',''+str(number_atoms+10))
        

        if build_indices[0] < 10:
            fg_bond_to_build = fg_bond_to_build.replace('z',' '+str(build_indices[0]))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('z',''+str(build_indices[0]))

        if len(build_indices) > 1:
            if build_indices[1] < 10:
                fg_bond_to_build = fg_bond_to_build.replace('y',' '+str(build_indices[1]))
            else:
                fg_bond_to_build = fg_bond_to_build.replace('y',''+str(build_indices[1]))


#        buildsite_to_letter_dict = {
#                                    0: 'u',
#                                    1: 'v',
#                                    2: 'w',
#                                    3: 'x',
#                                    4: 'y',
#                                    5: 'z',
#        }

#        print(build_indices)

#        for each_buildsite in range(len(build_indices)):
#            letter_to_replace = buildsite_to_letter_dict[each_buildsite]
#            print(letter_to_replace)
#            print(build_indices[each_buildsite])
#            if build_indices[each_buildsite] < 10:
#                fg_bond_to_build = fg_bond_to_build.replace(letter_to_replace,' '+str(build_indices[each_buildsite]))
#            else:
#                fg_bond_to_build = fg_bond_to_build.replace(letter_to_replace,''+str(build_indices[each_buildsite]))
        

        new_mol_block = ''
        for each_line_index, each_line in enumerate(mol_block.split('\n')):
            
            if each_line_index == 3:
                number_bonds = int(each_line[4:6])
                new_number_atoms = number_atoms + self.ATTACH_GROUP_NUMB_ATOMSBONDS[0]
                new_number_bonds = number_bonds + self.ATTACH_GROUP_NUMB_ATOMSBONDS[1]


                if new_number_atoms > 9 and new_number_bonds > 9:
                    each_line = each_line.replace(each_line[0:7], ' ' + str(new_number_atoms) + ' ' + str(new_number_bonds) + ' ' + each_line[7:])
                elif new_number_atoms < 10 and new_number_bonds > 9:
                    each_line = each_line.replace(each_line[0:7], '  ' + str(new_number_atoms) + ' ' + str(new_number_bonds) + ' ' + each_line[7:])
                elif new_number_atoms > 9 and new_number_bonds < 10:
                    each_line = each_line.replace(each_line[0:7], ' ' + str(new_number_atoms) + '  ' + str(new_number_bonds) + ' ' + each_line[7:])
                elif new_number_atoms < 10 and new_number_bonds < 10:
                    each_line = each_line.replace(each_line[0:7], '  ' + str(new_number_atoms) + '  ' + str(new_number_bonds) + ' ' +  each_line[7:])

            if each_line_index == 4 + number_atoms:
                each_line = each_line.replace(each_line,self.ATTACH_GROUP+fg_bond_to_build+each_line)

            if each_line_index >= 4 + number_atoms:   
                
                if self.REMOVE_BOND_FROM_RXNSITE is not None:
                    if ' '+str(indices_removebond[0]+1)+' 'in each_line[0:7] and ' ' + str(indices_removebond[1]+1) + ' ' in each_line[0:7]:
                        each_line = each_line.replace(each_line,each_line[0:7]+ ' '+str(int(each_line[7:9])-1)+each_line[9:])
                
                if self.ADD_BOND_TO_RXNSITE is not None:
                    if ' '+str(indices_addbond[0]+1)+' ' in each_line[0:7] and ' '+str(indices_addbond[1]+1)+' ' in each_line[0:7]:
                        each_line = each_line.replace(each_line, each_line[0:7] +' '+ str(int(each_line[7:9])+1) + each_line[9:])

            new_mol_block = new_mol_block + each_line + '\n'


        new_mol = Chem.MolFromMolBlock(new_mol_block)
        new_mol = Chem.AddHs(new_mol)

        # Embed the molecule
        AllChem.EmbedMolecule(new_mol, AllChem.ETKDG())
        # Perform geometry optimization
        AllChem.MMFFOptimizeMolecule(new_mol)

        temp_outmol_filename = 'temp2out.mol'
        temp_outxyz_filename = 'temp2out.xyz'
        

        # Save optimized molecule to a mol file
        Chem.MolToMolFile(new_mol, temp_outmol_filename)

        utils.convert_xyz2mol(temp_outxyz_filename,temp_outmol_filename,BACKWARDS=True)

        #read in out xyz file 
        opt_xyz_file = open(temp_outxyz_filename,mode='r')
        opt_prod_xyz_string = ''
        for line in opt_xyz_file:
            opt_prod_xyz_string = opt_prod_xyz_string + line 

        return opt_prod_xyz_string

    def RemoveBuild_LGFG(self):
        '''
        This removes the atom-labels associated with an atomic environment 
        that needs to be removed to make the product

            QM9_BOOL                - boolean if QM9 is being used to load with QM9 method so as to not have to list all the properties that QM9 has,
                                      otherwise load with AtomsData method (inconsequential)
            data                    - loading the entire db file using either the QM9 method or AtomsData method from schnetpack
            number_rxn              - keeps track of number of products/reactants found
            this_mol_idx            - integer running through the MOLECULES_RANGE of the db we are interested in exploring for reactants
                                      to a possible reaction
            MOLECULES_RANGE         - this is the range of db to scan to find reactants, ex. [0,10000]
            at, props               - a molecule's properties extracted from QM9
                                        at - molecular formula
                                        props - is a dictionary of properties which can be accessed with keywords 
            atomic_numbers          - list of atomic numbers in the molecule, from props 
            numb_atoms              - number of atoms in the molecule 
            xyz_positions           - coordinates of atoms in the molecule, from props
            react_xyz_string        - the possible target reactant's xyz string
            RemoveBuild_LGFG_
            SingleBuildSite.
            Write_InitXYZString     - function that writes the possible reactant's xyz string from 
                                      molecule props and returns it
            atomlabels_rxncenters   - the list of all integer labels at depth N that contain the reaction center
                                      (a list of lists one for each reaction center)
            
                                      


            indices_fg_found        - this list holds all the places where we found the atoms associated
                                      with the FG to be removed. Using the atomlabels_remove[depth][1], integer label, 
                                      to find the exact indices of FG to be removed
            indices_thisatomenv_
            found                   - some atomic environments may have two in a functional group, ex. the H's in NH2
                                      this is the list that contains them before they are all appended into indices_fg_found_above
            does_it_exist           - it is important to check that all the atoms of the FG are found 
                                      and not just one, before deleting the entire FG. This list holds booleans 'yes' and 'no'
                                      to determine if all the atoms of a FG have been found before the reaction can GO
                                      if all atoms don't exist the molecule is passed through without being counted for the reaction
            reactcenter_indices     - finds the index of the molecule at which the build atom is found
                                      THIS IS BEFORE FILTERING THE FG from the depthlabel_build_column_this_mol
                                      TO ENSURE that the build_index here can be used to extract the correct index and
                                      correct embeddings (if we do it after, then there will likely be a mismatch)
            all_reactcenter_indices - this is the list (from control) that will hold all the build_indices for extractembeddings
            filtered_build_
            depthlabel_this_mol     - filtered molecule labels of the building-depth with FG indices removed
            prodcenter_indices      - this list holds all the places where we can build the fg
                                      using the ATOMLABELS_BUILD[depth][1], integer label, to find the exact
                                      indices of the atoms that can build the fg off of (those that will lose an FG bonded to them)
                                      NOTE this has to be obtained AFTER finding indices of the functional group
                                      AND after filtering relevant_build_depth_labels_this_mol from the indices to delete
                                      if not done, the build indices may likely start ponting to the wrong atom,
                                      especially if the atoms deleted are before the build index then there is a frameshift
            all_prodcenter_indices  - the list that will hold all the build indices after FG has been removed
                                      important for embedding extraction of the target atom 
                                      FOR NOW THE TARGET ATOM FOR EMBEDDING IS THE ONE UNDERGOING THE MOST CHANGE 
                                      (WHERE THE LEAVING GROUP AND ATTACKING GROUP ARE HAPPENING)
                                      HOWEVER, IT IS POSSIBLE TO EXTRACT ANY OTHER TARGET ATOM BY INTRODUCING "all_targetembedding_indices"
                                      and look for the indices that have the label of the atom you want to extract (an atom in the FG!)
            removelg_xyz_string      - if all the atoms of the FG exist then we create the removelg_xyz_string
            RemoveBuild_LGFG_
            SingleBuildSite.
            RemoveAtoms_LG          - calling the function that removes atoms from react_xyz_string 
                                      based on indices found and returns removelg_xyz_string
            prod_xyz_string         - the product string which results from taking the FG and bonding to the 
                                      build_indices
            RemoveBuild_LGFG_
            SingleBuildSite.
            Build_FG                - calling the function that will bind the FG to the build_indices
                                      and return the prod_xyz_string      
            utils.generate_
            dbfromxyz               - utility function that converts xyz to db file so it is ready for embedding
                                      extraction                           
'''         

        #load the dataset
        #if the data is QM9 use schnetpack's QM9 loader as it is code-efficient
        #otherwise use AtomsData and make sure you define the available properties encoded in the db file
        if self.QM9_BOOL == True:
            data = QM9(self.DATASET_FILEPEATH,download=False,remove_uncharacterized=True)
        else:
            data = AtomsData(self.DATASET_FILEPEATH,available_properties=self.AVAILABLE_PROPERTIES)

        #initialize reaction counter
        number_rxn = 0

        #run through the dataset
        for this_mol_idx in range(self.MOLECULES_RANGE[0],self.MOLECULES_RANGE[1]):
            if this_mol_idx not in self.SKIP:
                #Some molecule will fail to optimize the product, you can skip them through the SKIP constant

                #load propeties of the molecule
                at, props = data.get_properties(this_mol_idx)
                atomic_numbers = props['_atomic_numbers'].detach().numpy()
                numb_atoms  = len(atomic_numbers)
                xyz_positions = props['_positions'].detach().numpy()
                
                #write the react_xyz
                self.react_xyz_string = RemoveBuild_LGFG_SingleBuildSite.Write_InitXYZString(self,atomic_numbers,numb_atoms,xyz_positions)

                #the atom-environment labels (at depth 5, or any sufficient depth) for each atom in this molecule
                envlabels_of_molecule = self.labels[self.labels['mol_idx'+str(self.DEPTH)] == this_mol_idx]['labeldepth'+str(self.DEPTH)].values

                #running the tool that finds the centers of reaction, IF THEY EXIST, their matched depth-structure labels, and their indices in the molecule
                does_eachrc_exist, rc_matched_labels, rc_indices = utils.find_if_reactioncenter_exists_exactconfig(self.RXN_CENTERS,self.atomlabels_rxncenters,envlabels_of_molecule)



                #IF any reaction center does not exist (or have the right count), or any leaving center does not exist (or have the right count)
                #we skip the molecule as a reactant for the fgtransform
                if 'no' in does_eachrc_exist:
                    #skip the molecule
                    pass
                #IF ALL reaction centers exist with the appropriate count, and all the leaving centers exist with the appropriate count, then this molecule
                #is a reactant! 
                else:
                    #reactant found, printing the molecule index
                    print('reactant found:',this_mol_idx)
                    
                    #does the environment exits?

                    if self.LEAVING_ATOMS is not None:
                        does_lc_exist, lc_matched_labels, lc_indices = utils.find_leavingatoms(self.LEAVING_ATOMS,self.atomlabels_leavinggroup,envlabels_of_molecule)

                    #if there is not leaving center you will need to bypass the does_each_lc_exist
                    #just to bypass the check that there is the correct leaving centers in this molecule
                    elif self.LEAVING_ATOMS == None:
                        does_lc_exist = ['yes']
                        lc_indices = []

                    if 'no' not in does_lc_exist:

                        #append to the list of all react indices so that we can extract embeddings of those atoms at the reaction center
                        self.all_reactcenterindices.append(rc_indices)
                        
                        #delete atoms of the molecule that correspond to the found leaving centers
                        envlabels_of_molecule_filteredlc = np.delete(envlabels_of_molecule,lc_indices)

                        #reassessing where matched react center labels are now that the leaving centers have been deleted in
                        #envlabels_of_molecule_filteredlc
                        pc_indices = np.where(np.isin(envlabels_of_molecule_filteredlc, rc_matched_labels))[0].tolist()      
                        #append to the list of all prod indices so that we can extract embeddings of those atoms at the prod centers
                        self.all_prodcenterindices.append(pc_indices)

                        #need to have these initialized incase one or other or both remains empty, which is always the case one of these
                        indices_removebond = []
                        indices_addbond = []
                        if self.REMOVE_BOND_FROM_RXNSITE is not None:
                            indices_removebond = utils.find_bondorderchange_atoms(self.REMOVE_BOND_FROM_RXNSITE,self.atomlabels_removebond,envlabels_of_molecule)
                        if self.ADD_BOND_TO_RXNSITE is not None:
                            indices_addbond = utils.find_bondorderchange_atoms(self.ADD_BOND_TO_RXNSITE,self.atomlabels_addbond,envlabels_of_molecule)

                        #remove the leaving group and give back the xyz string
                        #WITHOUT ANY OPTIMIZATION, this is a half-step to sometimes check atoms being removed without optimization
                        self.removelg_xyz_string = RemoveBuild_LGFG_SingleBuildSite.RemoveAtoms_LG(self,lc_indices)
                        #build the incoming functional group, add/remove any bonds at reaction site, the tool works with mol files, return the xyz string
                        self.prod_xyz_string = RemoveBuild_LGFG_SingleBuildSite.Build_FG(self,pc_indices,numb_atoms-len(lc_indices),indices_addbond,indices_removebond)
                        
                        #add number of reactions found
                        number_rxn = number_rxn + 1
                        
                        #write react, remove_lg, and prod into their respective files
                        self.react_dataset_xyzfile.write(self.react_xyz_string)
                        self.removelg_dataset_xyzfile.write(self.removelg_xyz_string)
                        self.prod_dataset_xyzfile.write(self.prod_xyz_string)

        #close the files of reactants, remove_lg, products
        self.react_dataset_xyzfile.close()
        self.removelg_dataset_xyzfile.close()  
        self.prod_dataset_xyzfile.close() 
        
        #write xyz into db files with available properties listed in input
        utils.generate_dbfromxyz(self.REACT_FILEPATH,self.AVAILABLE_PROPERTIES)
        utils.generate_dbfromxyz(self.REMOVELG_FILEPATH,self.AVAILABLE_PROPERTIES)
        utils.generate_dbfromxyz(self.PROD_FILEPATH,self.AVAILABLE_PROPERTIES)

        #return number of reactions, the react center indices, the prod center indices, which will be used to extract embeddings 
        #at the correct react center locations for reactant and products
        return number_rxn, self.all_reactcenterindices, self.all_prodcenterindices


                    
                