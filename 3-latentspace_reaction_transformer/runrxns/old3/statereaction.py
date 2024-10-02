from schnetpack.datasets import QM9
import numpy as np
from numpy import genfromtxt, savetxt

from fgtransform.utils import utils


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


class transform():
    '''
    date modified - 2024-02-16

    The main code that runs state reactions (A + B -> C + D, without modelling kinetics) on db files,
    it first scans for reactant-center labels and then performs the reaction, allows for optimization
    output are files of the reactants found and their corresponding products


        transform       - class that handles all types of transformations
                        through a variety of functions:
        
                        write_init_string   - writes the reactant file (init.xyz/db) once reactant label have been found
                        optimize            - runs MMFF94 optimization from RDKit on the product 
                        removeHfunction     - used by removeH, does the erasing from mol/xyz files removing H
                        removeH             - controls removal of a target labelled H
                        oxidalcsfunction    - used by oxidalcs, does the erasing from mol/xyz files thus oxidizing
                        oxidalcs            - control removal of H's associated with alcohols, their alpha carbon or their alpha carbon that is branched (C-CH-C)

    '''
    def __init__(self,dataset_filepath,targetlabelH_filepath,labelH1id,init_dataset_filepath,trans_dataset_filepath,n_molecules,available_properties):
        '''
        initializing the input variables that control the entire code,
        AND the build variables of the code which ultimately 
        construct the reactant init.xyz and product files trans.xyz

            control variables: 
                dataset_filepath            - the db filepath to scan for possible reactants with reactants label
                targetlabelH_filepath       - the label filepath associated with the db file
                labeH1id                    - the target H label that identifies the functional group for the reaction
                init_dataset_filepath       -
                trans_dataset_filepath      -
                n_molecules                 -
                available_properties        -

            build variables:   
                init_dataset_xyzfilepath    -
                trans_dataset_xyzfilepath   -
        '''

        #control variables
        self.dataset_filepath = dataset_filepath
        self.targetlabelH_filepath = targetlabelH_filepath
        self.labelH1id = labelH1id
        self.init_dataset_xyzfilepath = init_dataset_filepath
        self.trans_dataset_xyzfilepath = trans_dataset_filepath
        self.n_molecules = n_molecules
        self.available_properties = available_properties

        #build variables
        self.init_dataset_xyzfile = open(self.init_dataset_xyzfilepath, mode='w')
        self.trans_dataset_xyzfile = open(self.trans_dataset_xyzfilepath, mode='w')


    def write_init_string(self,atomic_numbers,numb_atoms,xyz_positions):
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



    #to use RDKit to optimize either init or trans structure
    def optimize(self,trans_xyz_string):

        #write a temp xyz string
        temp_xyz_filename = 'temp2.xyz'
        temp_xyzfile = open(temp_xyz_filename,mode='w')
        temp_xyzfile.write(trans_xyz_string)
        temp_xyzfile.close()

        temp_mol_filename = 'temp2.mol'
        utils.convertxyz_to_mol(temp_xyz_filename,temp_mol_filename)

        #
        mol = Chem.MolFromMolFile(temp_mol_filename)
        
        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Embed the molecule
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        # Perform geometry optimization
        AllChem.MMFFOptimizeMolecule(mol)

        temp_outmol_filename = 'temp2out.mol'
        temp_outxyz_filename = 'temp2out.xyz'
        # Save optimized molecule to a mol file
        Chem.MolToMolFile(mol, temp_outmol_filename)

        utils.convertxyz_to_mol(temp_outxyz_filename,temp_outmol_filename,backwards=True)

        #read in out xyz file 
        opt_xyz_file = open(temp_outxyz_filename,mode='r')
        opt_trans_xyz_string = ''
        for line in opt_xyz_file:
            opt_trans_xyz_string = opt_trans_xyz_string + line 

        return opt_trans_xyz_string


    #writing the transformed file
    def removeHfunction(self,H1_index_in_xyzmol):
        countH_inxyz = 0

        trans_xyz_string = ''
        for this_line_index, this_line in enumerate(self.init_xyz_string.split('\n')):
            if this_line_index == 0:
                this_line = this_line.replace(this_line,'\n'+str(int(this_line)-1)+'\n')

            if 'H' in this_line:
                if countH_inxyz == H1_index_in_xyzmol[0]:
                    this_line = this_line.replace(this_line,'')
                countH_inxyz = countH_inxyz + 1
        
            trans_xyz_string = trans_xyz_string + this_line +'\n'
        
        #remove empty lines
        lines = trans_xyz_string.split('\n')
        trans_xyz_nonempty_lines = [line for line in lines if line.strip()]

        # Join the non-empty lines back into a string
        trans_xyz_string2 = '\n'.join(trans_xyz_nonempty_lines)+'\n'
        return trans_xyz_string2

    def oxidalcsfunction(self,H1_index_in_xyzmol,H2_index_in_xyzmol):
        countH_inxyz = 0

        trans_xyz_string = ''
        for this_line_index, this_line in enumerate(self.init_xyz_string.split('\n')):
            if this_line_index == 0:
                this_line = this_line.replace(this_line,'\n'+str(int(this_line)-2)+'\n')

            if 'H' in this_line:
                if countH_inxyz == H1_index_in_xyzmol[0]:
                    this_line = this_line.replace(this_line,'')
                if countH_inxyz == H2_index_in_xyzmol[0]:
                    this_line = this_line.replace(this_line,'')
                countH_inxyz = countH_inxyz + 1
        
            trans_xyz_string = trans_xyz_string + this_line +'\n'
        
        #remove empty lines
        lines = trans_xyz_string.split('\n')
        trans_xyz_nonempty_lines = [line for line in lines if line.strip()]

        # Join the non-empty lines back into a string
        trans_xyz_string2 = '\n'.join(trans_xyz_nonempty_lines)+'\n'
        return trans_xyz_string2

    def removeH(self):
        
        #generate qm9data and H label data for qm9
        qm9data = QM9(self.dataset_filepath,download=False,remove_uncharacterized=True)
        labelH = genfromtxt(self.targetlabelH_filepath,delimiter=' ',encoding='utf-8-sig')

        number_trans = 0
        for this_molecule in range(self.n_molecules[0],self.n_molecules[1]):

            #load molecule properties
            at, props = qm9data.get_properties(this_molecule)
            xyz_positions = props['_positions'].detach().numpy()
            atomic_numbers = props['_atomic_numbers'].detach().numpy()
            numb_atoms = len(props['_atomic_numbers'].detach().numpy())

#            H_index  = oxidation.indexer(props,each_molecule,labelH)

            #only if ONE oxygen exists in the molecule 
            #(FOR now, difficult to know which HO associates with which other H for oxidation unless its one)
            H1_found = False
            if np.count_nonzero(atomic_numbers==8) == 1: 

                #fine the first occurence H_index where the molecule is (contained in the label_file)
                H_index = np.where(labelH[:,1] == this_molecule)

                #Find out the two H_indices that contain the two H's 
                # undergoing oxidation
                for this_H_in_molecule in range(len(H_index[0])):
                    if labelH[H_index[0][this_H_in_molecule]][0] == self.labelH1id:


                        H1_index_in_xyzmol = np.where(H_index[0] == H_index[0][this_H_in_molecule])[0]
                        H1_found = True


                if H1_found == True:
                    number_trans = number_trans + 1
                    #make xyz string and write scratch xyz file
                    self.init_xyz_string = transform.write_init_string(self,atomic_numbers,numb_atoms,xyz_positions)

#                   scratch_xyz_filename = self.scratch_file + '.xyz' 
                    self.init_dataset_xyzfile.write(self.init_xyz_string)

                    #make transformed xyz
                    trans_xyz_string = transform.removeHfunction(self,H1_index_in_xyzmol)
                    self.trans_dataset_xyzfile.write(trans_xyz_string)

                else:
                    pass

        self.init_dataset_xyzfile.close()
        self.trans_dataset_xyzfile.close()
        

        utils.generate(self.init_dataset_xyzfilepath,self.available_properties,self.n_molecules[1])
        utils.generate(self.trans_dataset_xyzfilepath,self.available_properties,self.n_molecules[1])

        return number_trans
    

    def oxidalcs(self,targetlabelH_filepathacarb,labelH2id,labelH3id,optimize):

        #generate qm9data and H label data for qm9
        qm9data = QM9(self.dataset_filepath,download=False,remove_uncharacterized=True)
        labelhydroxH = genfromtxt(self.targetlabelH_filepath,delimiter=' ',encoding='utf-8-sig')
        labelacarbH = genfromtxt(targetlabelH_filepathacarb,delimiter=' ',encoding='utf-8-sig')

        number_trans = 0


        for this_molecule in range(self.n_molecules[0],self.n_molecules[1]):

            #loadbar
            if this_molecule % 1000 == 0:
                print(this_molecule)

            #load molecule properties
            at, props = qm9data.get_properties(this_molecule)
            xyz_positions = props['_positions'].detach().numpy()
            atomic_numbers = props['_atomic_numbers'].detach().numpy()
            numb_atoms = len(props['_atomic_numbers'].detach().numpy())

            H1_found = False
            H2_found = False
            if np.count_nonzero(atomic_numbers==8) == 1: 

                #fine the first occurence H_index where the molecule is (contained in the label_file)
                H_indexhydroxy = np.where(labelhydroxH[:,1] == this_molecule)
                H_indexacarb = np.where(labelacarbH[:,1] == this_molecule) 

                #Find out the two H_indices that contain the two H's 
                # undergoing oxidation
                for this_H_in_molecule in range(len(H_indexhydroxy[0])):
                    if labelhydroxH[H_indexhydroxy[0][this_H_in_molecule]][0] == self.labelH1id:

                        H1_index_in_xyzmol = np.where(H_indexhydroxy[0] == H_indexhydroxy[0][this_H_in_molecule])[0]
                        H1_found = True
                    
                    if labelacarbH[H_indexacarb[0][this_H_in_molecule]][0] == labelH2id or labelacarbH[H_indexacarb[0][this_H_in_molecule]][0] == labelH3id:
                        H2_index_in_xyzmol = np.where(H_indexacarb[0] == H_indexacarb[0][this_H_in_molecule])[0]
                        H2_found = True

                if H1_found == True and H2_found == True:
                    trans_found = False
                    try:
                        self.init_xyz_string = transform.write_init_string(self,atomic_numbers,numb_atoms,xyz_positions)
                        #make transformed xyz
                        trans_xyz_string = transform.oxidalcsfunction(self,H1_index_in_xyzmol,H2_index_in_xyzmol)
                        if optimize == True: 
                            trans_xyz_string = transform.optimize(self,trans_xyz_string)
                            self.trans_dataset_xyzfile.write(trans_xyz_string)
                            trans_found = True
                        elif optimize == False:
                            self.trans_dataset_xyzfile.write(trans_xyz_string)
                            trans_found = True

                    except ZeroDivisionError:
                        print(f"Error at index {this_molecule}: Division by zero")
                    except IndexError:
                        print(f"Error at index {this_molecule}: Index out of bounds")
                    except ValueError:
                        print(f"Error at index {this_molecule}: Value error")
                    except Exception as e:
                        print(f"Error at index {this_molecule}: {e}")

                    if trans_found == True:
                        number_trans = number_trans + 1
                        #make xyz string and write scratch xyz file
    #                   scratch_xyz_filename = self.scratch_file + '.xyz' 
                        self.init_dataset_xyzfile.write(self.init_xyz_string)



        self.init_dataset_xyzfile.close()
        self.trans_dataset_xyzfile.close()

        utils.generate(self.init_dataset_xyzfilepath,self.available_properties,self.n_molecules[1])
        utils.generate(self.trans_dataset_xyzfilepath,self.available_properties,self.n_molecules[1])

        return number_trans
    