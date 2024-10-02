
from schnetpack.datasets import QM9
from schnetpack import AtomsData
import pandas as pd
import csv
import numpy as np

from fgtransform.utils import utils

from rdkit import Chem
from rdkit.Chem import AllChem

import warnings

warnings.filterwarnings("ignore")

'''
Last Updated: 2024-03-01

'''

class RemoveBuild_LGFG_SingleBuildSite():
    '''
    transforms reactants to products by cutting/pasting depth-labelled atoms.
    FOR NOW the code ONLY handles one build site. Having two build sites might be problematic
    in terms of symmetry. 
    
        functions in the class

            __init__             - initializes the control and build variables
            Write_InitXYZString  - write the reactant if atomic environment is found in the db file
            RemoveAtoms_LG       - function that removes atoms of the atomic environment from an xyz file
                                   returns the resulting xyz without any optimization for analysis of embeddings
                                   of this intermediatte step (before adding the functional group)
            Build_FG             - function that build the functional group in the atomic environment 
                                   and optimizes it using RDKit's MMFF94, 
                                   NOTE requires conversion to mol file with obabel because functional groups
                                   have to be tacked on to the correct place (build_index) 
            RemoveBuild_LGFG     - the function that works other functions above to find and remove atom labels 
                                   associated with an atomic environment and then build a functional group
                                   instead of the removed functional group
    '''

    def __init__(self, DATASET_FILEPEATH,LABELS_FILEDIR,ATOMLABELS_REMOVE,ATOMLABELS_BUILD,INIT_DATASET_FILEPATH,REMOVE_DATASET_FILEPATH,TRANS_DATASET_FILEPATH,MOLECULES_RANGE,AVAILABLE_PROPERTIES,QM9,FG_TO_BUILD,FG_BOND_TO_BUILD,ADD_ATOMS_BONDS,N_DEPTHS,ONLY_ONE,FG_ADD_DOUBLE_BOND,ADD_DOUBLE_BOND_BW,SKIP):
        '''
        Initializes the control  & build variables:

            control:

                DATASET_FILEPEATH           - the db file that has been labelled
                QM9                         - a boolean depending on qm9.db is being used because to load all the properties of that one
                                              we use the QM9 method otherwise the AtomsData method from schnetpack
                LABELS_FILEDIR              - in a specified directory, you should have all the labels at every depth, "0.csv", "1.csv"
                                              all of these will be stacked horizontally and accessed with pandas 
                                              this will allow ease of access to labels at all depths so
                                              that the algorithm can pick out the correct atom to remove using their environment label and
                                              depth specified
                ATOMLABELS_REMOVE           - a 2-column data matrix that defines the FULL ATOMIC ENVIRONMENT! 

                                              depth1  |  label_to_remove1
                                              depth2  |  label_to_remove2
                                               
                                              ex.
                                              [['labeldpeth3',31],['labeldepth2',24]]

                                              this helps the code identify which depth the environment label to look for
                                              and the environment label itself, in order to remove the atom in the correct environment

                                              often an environment is multiple atoms at different depths to specify it fully
                                              without being too restrictive, for instance, removin H-O-C- requires H depth level 2
                                              to ensure that there is carbon, but O depth level 1 so that we are not too restrictive
                                              on what is after the carbon. All these environments should be pre-labelled (using the autolabel code) 
                                              and can be found in LABELS_FILEDIR


                ATOMLABELS_BUILD            - a list that tells the algorithm which atomlabel to build the functional group on
                                              this requires also giving the depth of the atomlabel and the label of the build atom
                                              ex.

                                              ['labeldepth1',24]
                
                MOLECULES_RANGE             - range of molecules from the db file to scan for atomic environments and thus allow 
                                              through the algorithm for transformation (removing atomic environments, adding functional groups)
                                              
                                              NOTE this range needs to be labeled by autolabel in the label_filedir!
                
                AVAILABLE_PROPERTIES        - this is required for making the db files from xyz files,
                                              db files comes with associated properties, 
                                              for now this is a filler! and can be anything! just keep note of what name it is so the db can be accessed later
                FG_TO_BUILD                 - this is mol string information about the functional group to build
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
                ADD_ATOMS_BONDS             - this specifies to the algorithm how many atoms and bonds this functional group adds 
                                              so it can adjust the mol file accordingly

                                              ex.
                                              [2,2]    - add two atoms, add two bonds
                FG_BOND_TO_BUILD            - a template for the bonding section of the mol file that will replace alphabetical symbols
                                              with the correct connections
                                            
                                              ex.
                                              fg_OH_bond = ' a f  1  0\n b a  1  0\n'

                                              the code knows that 'a' means index 1, 'b' means 2, and so on
                                              with the exception of 'f' which is replaced with the bond_index of the functional group
                                              found by using the atombuild_label
                                              this is kept as a template copy and the Build_FG method copies it 
                                              replace the alphabets with the correct connecting indices

                N_DEPTHS                    - the number of depths in the label_filedir
                ONLY_ONE                    - some removals require only one of the atomic environment to be removed not both (i.e CH3CH2)
                                              some removals require all atomic environments to be removed (C=O-(NH2))
                                              this boolean controls this 
                FG_ADD_DOUBLE_BOND          - boolean to decide to add double bond as a way to esnure a double bond is created
                ADD_DOUBLE_BOND_BW          - the labels of the atomic environments to create the double bond between               
                SKIP                        - list of indices to skip because they failed at geometry optimization
                all_label_depths            - this is a horizontal concatenation of all the label depths from the
                                              label_filedir, allows ease of access to all the depthlabels at once
                                              so that when 'label1' is specified that column for depth1 is picked out
                


            build:

                INIT_DATASET_FILEPATH       - the filepath to create the reactants file
                REMOVE_DATASET_FILEPATH     - the filepath to create the removal of leaving group ONLY 
                                              without optimization
                TRANS_DATASET_FILEPATH      - the filepath to create the products file 
                                              with optimization
                init_dataset_xyzfile        - open the init file in write mode
                remove_dataset_xyzfile      - open the rem file in write mode
                trans_dataset_xyzfile       - open the trans file in write mode
                all_build_indices           - this is important to tag the build_index in a list so that
                                              the extractembeddings file can extract embeddings of
                                              the build atom 
                                              (where the leaving group is detaching and attacking group attacking, for example)
                all_build_initindices       - after removal of the leaving, the tagged build index may be different
                                              than before the leaving group, this ensures a tag before the leaving group
                                              for the embeddings of the reactant atom to be extracted

        '''

        #control variables
        self.DATASET_FILEPEATH = DATASET_FILEPEATH
        self.QM9 = QM9
        self.LABELS_FILEDIR = LABELS_FILEDIR
        self.ATOMLABELS_REMOVE  = ATOMLABELS_REMOVE 
        self.ATOMLABELS_BUILD   = ATOMLABELS_BUILD
        self.MOLECULES_RANGE = MOLECULES_RANGE
        self.AVAILABLE_PROPERTIES = AVAILABLE_PROPERTIES
        self.FG_TO_BUILD = FG_TO_BUILD
        self.ADD_ATOMS_BONDS = ADD_ATOMS_BONDS
        self.FG_BOND_TO_BUILD = FG_BOND_TO_BUILD
        self.N_DEPTHS = N_DEPTHS
        self.only_one = ONLY_ONE
        self.FG_ADD_DOUBLE_BOND = FG_ADD_DOUBLE_BOND
        self.ADD_DOUBLE_BOND_BW = ADD_DOUBLE_BOND_BW
        self.SKIP = SKIP

        #load and concatenate (horizontally) the labels at all depths from the autolabel filedirectory
        self.labels_alldepths = pd.DataFrame()
        for i in range(0,self.N_DEPTHS):
            file_path = f'{self.LABELS_FILEDIR}/{i}.csv'
            with open(file_path, mode='r') as file:
                csv = pd.read_csv(file,index_col=None)
                self.labels_alldepths = pd.concat([self.labels_alldepths, csv], axis=1)



        #build variables
        #these filepaths will be the same place where xyz will be converted to the db file
        self.init_dataset_xyzfilepath = INIT_DATASET_FILEPATH
        self.remove_dataset_xyzfilepath = REMOVE_DATASET_FILEPATH
        self.trans_dataset_xyzfilepath = TRANS_DATASET_FILEPATH
        self.init_dataset_xyzfile = open(self.init_dataset_xyzfilepath, mode='w')
        self.remove_dataset_xyzfile = open(self.remove_dataset_xyzfilepath,mode='w')
        self.trans_dataset_xyzfile = open(self.trans_dataset_xyzfilepath, mode='w')

        #save all_build_indices for each molecule,
        #there will be ONLY 1 per molecule, 
        #this is so that we can use it to find 
        #the embedding of the target build atom
        #for init_ it has to be found by label,
        # as build_label may be different before/after removal
        self.all_build_indices = []
        self.all_build_initindices = []
        
        #You can set your own target, for embedding extraction
        #it doesn't have to be the build_index of each molecule, for later
        #NEXT UPDATE
        self.all_targetextraction_initindices = []
        self.all_targetextraction_indices = []



    def Write_InitXYZString(self,atomic_numbers,numb_atoms,xyz_positions):
        '''
        Writes the reactant if the atomic environment has been found in the db file

            init_xyz_string         - the xyz string that will be written based on xyz
                                      information about the molecules
            numb_atom               - number of atoms, info comes from the db file
            xyz_positions           - xyz coordinates, info comes from the db file
            each_atom               - an integer (0 -> N_atoms-1) for each atom in molecule
            atomic_numbers          - atomic numbers list, info comes from the db file
        ''' 
        
        init_xyz_string = ''
        init_xyz_string = init_xyz_string + str(numb_atoms) + '\n0.0000 \n'
        for each_atom in range(numb_atoms):
            if atomic_numbers[each_atom] == 1:
                init_xyz_string = init_xyz_string + 'H' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 6:
                init_xyz_string = init_xyz_string + 'C' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 7:
                init_xyz_string = init_xyz_string + 'N' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 8:
                init_xyz_string = init_xyz_string + 'O' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
            if atomic_numbers[each_atom] == 9:
                init_xyz_string = init_xyz_string + 'F' + '  ' + str(xyz_positions[each_atom][0]) + '  ' + str(xyz_positions[each_atom][1]) + '  ' + str(xyz_positions[each_atom][2]) + '\n'
        return init_xyz_string

    def RemoveAtoms_LG(self,indices_remove):
        '''
        function that removes the leaving group from the molecule 
        (as specified by ATOMLABELS_REMOVE)

            remove_xyz_string           - the resulting xyz_string which will have the leaving group removed
            count_xyzline               - initial counter for how many xyz lines have been passed
            this_line_index             - indices of lines in the xyz file
            this_line                   - lines in the xyz file
            init_xyz_string             - the initial xyz string which needs to be modifed by this function 
                                          (by removing leaving group)
            indices_remove              - the indices of the atoms to remove 
                                          (obtained previously from using the ATOMLABELS_REMOVE found for this molecule)
            remove_xyz_withemptylines   - the lines of the removed_xyz_string which may have
                                          some lines empty, now that the leaving group has been removed
            remove_xyz_stringnonempty_
            lines_list                  - a list that has all the lines in a 
            remove_xyz_string_noempty   - it is important for this to work 
                                          that the resulting string is not copied on top
                                          of remove_xyz_string, we use this copy to join the lines_list 
                                          above with '\n' thus giving us the xyz with LG removed
                                          WITH NO EMPTY LINES

                                          IT WILL NOT WORK IF YOU USE ORIGINAL BECAUSE
                                          IT WILL THINK TO JOIN '\n' the lines from the the original as well
        '''

        remove_xyz_string = ''
        count_xyzline = 0

        for this_line_index, this_line in enumerate(self.init_xyz_string.split('\n')):
            
            if this_line_index == 0:
                this_line = this_line.replace(this_line,'\n'+str(int(this_line)-len(indices_remove)) + '\n')
            
            if this_line_index > 1: 
                if count_xyzline in indices_remove:
                    this_line = this_line.replace(this_line,'')
                count_xyzline = count_xyzline + 1

            remove_xyz_string = remove_xyz_string + this_line + '\n'

        #remove empty lines
        remove_xyz_withemptylines = remove_xyz_string.split('\n')
        remove_xyz_stringnonempty_lines = [line for line in remove_xyz_withemptylines if line.strip()]
        # Join the non-empty lines back into a string
        remove_xyz_string_noempty = '\n'.join(remove_xyz_stringnonempty_lines)+'\n'

        return remove_xyz_string_noempty

    def Build_FG(self,build_index,number_atoms,indices_double):
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
            opt_trans_xyz_string        - writing the optimized final structure as a string and returning that
        '''

        build_index  = build_index[0]+1 #to match the indexing of mol file bonding section which starts at 1

        temp_xyz_filename = 'temp.xyz'
        temp_xyz_file = open(temp_xyz_filename, mode='w')
        temp_xyz_file.write(self.remove_xyz_string)
        temp_xyz_file.close()

        temp_mol_filename = 'temp.mol'
        utils.convert_xyz2mol(temp_xyz_filename,temp_mol_filename)

        mol = Chem.MolFromMolFile(temp_mol_filename,sanitize=False,removeHs=False)
#        mol = Chem.AddHs(mol)

        mol_block = Chem.MolToMolBlock(mol)
    
        fg_bond_to_build = self.FG_BOND_TO_BUILD 

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
        
        if build_index < 10:
            fg_bond_to_build = fg_bond_to_build.replace('z',' '+str(build_index))
        else:
            fg_bond_to_build = fg_bond_to_build.replace('z',''+str(build_index))

        new_mol_block = ''
        for each_line_index, each_line in enumerate(mol_block.split('\n')):
            
            if each_line_index == 3:
                number_bonds = int(each_line[4:6])
                new_number_atoms = number_atoms + self.ADD_ATOMS_BONDS[0]
                new_number_bonds = number_bonds + self.ADD_ATOMS_BONDS[1]


                if new_number_atoms > 9 and new_number_bonds > 9:
                    each_line = each_line.replace(each_line[0:7], ' ' + str(new_number_atoms) + ' ' + str(new_number_bonds) + ' ' + each_line[7:])
                elif new_number_atoms < 10 and new_number_bonds > 9:
                    each_line = each_line.replace(each_line[0:7], '  ' + str(new_number_atoms) + ' ' + str(new_number_bonds) + ' ' + each_line[7:])
                elif new_number_atoms > 9 and new_number_bonds < 10:
                    each_line = each_line.replace(each_line[0:7], ' ' + str(new_number_atoms) + '  ' + str(new_number_bonds) + ' ' + each_line[7:])
                elif new_number_atoms < 10 and new_number_bonds < 10:
                    each_line = each_line.replace(each_line[0:7], '  ' + str(new_number_atoms) + '  ' + str(new_number_bonds) + ' ' +  each_line[7:])

            if each_line_index == 4 + number_atoms:
                each_line = each_line.replace(each_line,self.FG_TO_BUILD+fg_bond_to_build+each_line)

            if each_line_index >= 4 + number_atoms:   
                if self.FG_ADD_DOUBLE_BOND == True:
                    if ' '+str(indices_double[0])+' ' in each_line[0:7] and ' '+str(indices_double[1])+' ' in each_line:
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
        opt_trans_xyz_string = ''
        for line in opt_xyz_file:
            opt_trans_xyz_string = opt_trans_xyz_string + line 

        return opt_trans_xyz_string

    def RemoveBuild_LGFG(self):
        '''
        This removes the atom-labels associated with an atomic environment 
        that needs to be removed to make the product

            QM9                     - boolean if QM9 is being used to load with QM9 method so as to not have to list all the properties that QM9 has,
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
            init_xyz_string         - the possible target reactant's xyz string
            RemoveBuild_LGFG_
            SingleBuildSite.
            Write_InitXYZString     - function that writes the possible reactant's xyz string from 
                                      molecule props and returns it
            ATOMLABELS_REMOVE       - the labels of the fg to remove, might require multiple depths,
                                      MUST BE organized from largest depth to smallest 
                                      format: [['labeldepth3',44],['labeldepth2',33]]           
            remove_depthlabels_
            columns_this_mol        - Because each local functional group that needs to be removed might have multiple atoms
                                      each atom will likely require a different depth depending on its positioning within 
                                      the functional group of the reaction. 
                                      
                                      This matrix filters down labels_alldepths (from input control) 
                                      down to on only the depths required to pick out the atoms of the FG
                                      
                                      and also filters to only those labels of this_mol_idx (this molecule)
                                      
                                      (see input control for labels_alldepths and how that was concatenated from labeldepths file directory)
                                            
                                        self.labels_alldepths[self.labels_alldepths['mol_idx0'] == this_mol_idx][[self.ATOMLABELS_REMOVE[each_chosen_depth][0] for each_chosen_depth in range(len(self.ATOMLABELS_REMOVE))]]
        
                                                                ^^^ this gets the all the labels of this molecule,                       ^^^^ list comprehension to only load those depths (each_chosen_depth) involved in  
                                                                    NOTE it does not matter if you choose 'mol_idx0'                          the labeling of environment that needs to be removed, see ATOMLABELS_REMOVE     
                                                                    column or any other ('mol_idx1' .. where 0 and 1...                       [['labeldepth3',63]] as an example, 
                                                                    correspond  to depths), bceause there is a label for 
                                                                    each depth so all of them will be isolated for this molecule
                                                                    whichever mol_idx column is chosen. It is this way because
                                                                    all depth labels got horizontally stacked.
            ATOMLABELS_BUILD        - the depth/label of the build atom (where to bond the fg), ex. ['labeldepth1',label#]
            build_depthlabel_
            column_this_mol         - This is the labels of the molecule using the depth/label
                                      that correctly identifies that target build index
                                      it tells the algorithm to bond the functional group to the 
                                      correct location. The target atom has only one depth 
                                      to select out from all the depth, this selects that depth out
                                      which is specified by ATOMLABELS_BUILD
                                      NOTE that this will be filtered later by deleting all those indices of molecule 
                                      involved in the FG that was removed, this is to allow us to find the build indices 
                                      after the removal has taken place (otherwise there will likely be a mismatch)
                                      NOTE ALSO this is used to determine if there is 
                                      EXACTLY 1 of the bond index FOUND in molecule, 
                                      ONLY those molecules are allowed to react,
                                      the issue with dealing with symmetry of having two places to bind is more difficult for now
            indices_fg_found        - this list holds all the places where we found the atoms associated
                                      with the FG to be removed. Using the ATOMLABELS_REMOVE[depth][1], integer label, 
                                      to find the exact indices of FG to be removed
            only_one                - boolean to control when to remove all atomic environments of a type found and when to remove just one
                                      some cases of LG, multiple environments of the same type might have to be removed (i.e NH2 in amide to acids, the H2)
                                      some cases of LG, onle one of environement of the same type might have to be removed
                                      (i.e CH3CH2 oxidation, need to only remove one H on each)
            indices_thisatomenv_
            found                   - some atomic environments may have two in a functional group, ex. the H's in NH2
                                      this is the list that contains them before they are all appended into indices_fg_found_above
            does_it_exist           - it is important to check that all the atoms of the FG are found 
                                      and not just one, before deleting the entire FG. This list holds booleans 'yes' and 'no'
                                      to determine if all the atoms of a FG have been found before the reaction can GO
                                      if all atoms don't exist the molecule is passed through without being counted for the reaction
            initindex_build_found   - finds the index of the molecule at which the build atom is found
                                      THIS IS BEFORE FILTERING THE FG from the depthlabel_build_column_this_mol
                                      TO ENSURE that the build_index here can be used to extract the correct index and
                                      correct embeddings (if we do it after, then there will likely be a mismatch)
            all_build_initindices   - this is the list (from control) that will hold all the build_indices for extractembeddings
            filtered_build_
            depthlabel_this_mol     - filtered molecule labels of the building-depth with FG indices removed
            index_build             - this list holds all the places where we can build the fg
                                      using the ATOMLABELS_BUILD[depth][1], integer label, to find the exact
                                      indices of the atoms that can build the fg off of (those that will lose an FG bonded to them)
                                      NOTE this has to be obtained AFTER finding indices of the functional group
                                      AND after filtering relevant_build_depth_labels_this_mol from the indices to delete
                                      if not done, the build indices may likely start ponting to the wrong atom,
                                      especially if the atoms deleted are before the build index then there is a frameshift
            all_build_indices       - the list that will hold all the build indices after FG has been removed
                                      important for embedding extraction of the target atom 
                                      FOR NOW THE TARGET ATOM FOR EMBEDDING IS THE ONE UNDERGOING THE MOST CHANGE 
                                      (WHERE THE LEAVING GROUP AND ATTACKING GROUP ARE HAPPENING)
                                      HOWEVER, IT IS POSSIBLE TO EXTRACT ANY OTHER TARGET ATOM BY INTRODUCING "all_targetembedding_indices"
                                      and look for the indices that have the label of the atom you want to extract (an atom in the FG!)
            remove_xyz_string       - if all the atoms of the FG exist then we create the remove_xyz_string
            RemoveBuild_LGFG_
            SingleBuildSite.
            RemoveAtoms_LG          - calling the function that removes atoms from init_xyz_string 
                                      based on indices found and returns remove_xyz_string
            trans_xyz_string        - the product string which results from taking the FG and bonding to the 
                                      build_indices
            RemoveBuild_LGFG_
            SingleBuildSite.
            Build_FG                - calling the function that will bind the FG to the build_indices
                                      and return the trans_xyz_string      
            utils.generate_
            dbfromxyz               - utility function that converts xyz to db file so it is ready for embedding
                                      extraction                           
'''         

        #load the dataset
        if self.QM9 == True:
            data = QM9(self.DATASET_FILEPEATH,download=False,remove_uncharacterized=True)
        else:
            data = AtomsData(self.DATASET_FILEPEATH,AVAILABLE_PROPERTIES=self.AVAILABLE_PROPERTIES)

        #initialize reaction counter
        number_rxn = 0

        #run through the dataset
        for this_mol_idx in range(self.MOLECULES_RANGE[0],self.MOLECULES_RANGE[1]):
            if this_mol_idx not in self.SKIP:
                at, props = data.get_properties(this_mol_idx)
                atomic_numbers = props['_atomic_numbers'].detach().numpy()
                numb_atoms  = len(atomic_numbers)
                xyz_positions = props['_positions'].detach().numpy()
                #write the init_xyz
                self.init_xyz_string = RemoveBuild_LGFG_SingleBuildSite.Write_InitXYZString(self,atomic_numbers,numb_atoms,xyz_positions)

    #           #columns for removal
                remove_depthlabels_columns_this_mol = self.labels_alldepths[self.labels_alldepths['mol_idx0'] == this_mol_idx][[self.ATOMLABELS_REMOVE[each_chosen_depth][0] for each_chosen_depth in range(len(self.ATOMLABELS_REMOVE))]]
                #columns for building
                build_depthlabel_column_this_mol = self.labels_alldepths[self.labels_alldepths['mol_idx0']==this_mol_idx][self.ATOMLABELS_BUILD[0]]

                #does the environment exits?
                does_it_exist = []
                indices_fg_found = []
                for each_depth in range(len(self.ATOMLABELS_REMOVE)):

                    if self.ATOMLABELS_REMOVE[each_depth][1] in remove_depthlabels_columns_this_mol.iloc[:,each_depth].tolist():
                        indices_thisatomenv_found = np.where(remove_depthlabels_columns_this_mol.iloc[:,each_depth].values == self.ATOMLABELS_REMOVE[each_depth][1])[0]

                        for each_idx in range(len(indices_thisatomenv_found)):

                            indices_fg_found.append(indices_thisatomenv_found[each_idx])
                            if self.only_one == True: 
                                break

                        does_it_exist.append('yes')
                        
                    else:
                        does_it_exist.append('no')

                if 'no' in does_it_exist:
                    pass
                else:
                    #
                    matched_label = None
                    for label in self.ATOMLABELS_BUILD[1:]:  # Start from index 1
#                        print(label)
                        if label in build_depthlabel_column_this_mol.tolist():
                            if matched_label is not None:  # Another match found, so not unique
                                matched_label = None
                                break
                            else:
                                matched_label = label

                    if np.count_nonzero(build_depthlabel_column_this_mol.values[build_depthlabel_column_this_mol == matched_label]) == 1:

#                        print(matched_label)

                        print(this_mol_idx)

                        #THIS HAS TO BE DONE HERE!!! 
                        #now you can find the build indices using the filtered labels
                        initindex_build_found = np.where(build_depthlabel_column_this_mol == matched_label)[0]
                        self.all_build_initindices.append(initindex_build_found.tolist())
                        
                        #THIS HAS TO BE DONE HERE!!!! After you find indices, 
                        #need to delete all indices
                        filtered_build_depthlabel_this_mol = np.delete(build_depthlabel_column_this_mol.values,indices_fg_found,axis=0)

                        #THIS HAS TO BE DONE HERE!!! 
                        #now you can find the build indices using the filtered labels
                        index_build = (np.where(filtered_build_depthlabel_this_mol == matched_label)[0])
                        #save build indices for embedding extraction
                        self.all_build_indices.append(index_build.tolist())

                        indices_double = []
                        if self.FG_ADD_DOUBLE_BOND == True:
                            indices_double = [np.where(filtered_build_depthlabel_this_mol == self.ADD_DOUBLE_BOND_BW[0])[0][0]+1, np.where(filtered_build_depthlabel_this_mol == self.ADD_DOUBLE_BOND_BW[1])[0][0]+1]

    #                    indices_build_foundin_init = np.where(relevant_remove_depthlabels_this_mol)

                        self.remove_xyz_string = RemoveBuild_LGFG_SingleBuildSite.RemoveAtoms_LG(self,indices_fg_found)
                            
                        #need to remove/pop all indices found from the relevant
                                
                        self.trans_xyz_string = RemoveBuild_LGFG_SingleBuildSite.Build_FG(self,index_build,numb_atoms-len(indices_fg_found),indices_double)
                        number_rxn = number_rxn + 1


                        self.init_dataset_xyzfile.write(self.init_xyz_string)
                        self.remove_dataset_xyzfile.write(self.remove_xyz_string)
                        self.trans_dataset_xyzfile.write(self.trans_xyz_string)

        self.init_dataset_xyzfile.close()
        self.remove_dataset_xyzfile.close()  
        self.trans_dataset_xyzfile.close() 
        utils.generate_dbfromxyz(self.init_dataset_xyzfilepath,self.AVAILABLE_PROPERTIES)
        utils.generate_dbfromxyz(self.remove_dataset_xyzfilepath,self.AVAILABLE_PROPERTIES)
        utils.generate_dbfromxyz(self.trans_dataset_xyzfilepath,self.AVAILABLE_PROPERTIES)

        return number_rxn, self.all_build_indices, self.all_build_initindices
