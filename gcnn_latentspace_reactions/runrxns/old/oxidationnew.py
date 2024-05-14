from schnetpack.datasets import QM9
import numpy as np
from numpy import genfromtxt, savetxt

from fgtransform.utils import utils

class oxidation():
    def __init__(self,dataset_filepath,targetlabelH_filepath,labelH1id,labelH2id,init_dataset_filepath,trans_dataset_filepath,n_molecules,atom_num,available_properties):
        self.dataset_filepath = dataset_filepath
        self.targetlabelH_filepath = targetlabelH_filepath
        self.labelH1id = labelH1id
        self.labelH2id = labelH2id
        self.init_dataset_xyzfilepath = init_dataset_filepath
        self.trans_dataset_xyzfilepath = trans_dataset_filepath
        self.n_molecules = n_molecules
        self.atom_num = atom_num

        self.init_dataset_xyzfile = open(self.init_dataset_xyzfilepath, mode='w')
        self.trans_dataset_xyzfile = open(self.trans_dataset_xyzfilepath, mode='w')

        self.available_properties = available_properties

    def noopt(self):
        
        #generate qm9data and H label data for qm9
        qm9data = QM9(self.dataset_filepath,download=False,remove_uncharacterized=True)
        labelH = genfromtxt(self.targetlabelH_filepath,delimiter=' ',encoding='utf-8-sig')

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
            H2_found = False
            if np.count_nonzero(atomic_numbers==8) == 1: 

                #fine the first occurence H_index where the molecule is (contained in the label_file)
                H_index = np.where(labelH[:,1] == this_molecule)

                #Find out the two H_indices that contain the two H's 
                # undergoing oxidation
                for this_H_in_molecule in range(len(H_index[0])):
                    if labelH[H_index[0][this_H_in_molecule]][0] == self.labelH1id:


                        H1_index_in_xyzmol = np.where(H_index[0] == H_index[0][this_H_in_molecule])[0]
                        H1_found = True

                    if labelH[H_index[0][this_H_in_molecule]][0] == self.labelH2id:
                    
                        H2_index_in_xyzmol = np.where(H_index[0] == H_index[0][this_H_in_molecule])[0]
                        H2_found = True


                if H1_found and H2_found == True:
                    print('molecule_idx',this_molecule)
                    #make xyz string and write scratch xyz file
                    init_xyz_string = ''
                    init_xyz_string = init_xyz_string + str(numb_atoms) + '\n0.0000 \n'
                    for each_atom in range(numb_atoms):
                        if atomic_numbers[each_atom] == 1:
                            init_xyz_string = init_xyz_string + 'H' + ' ' + str(xyz_positions[each_atom][0]) + ' ' + str(xyz_positions[each_atom][1]) + ' ' + str(xyz_positions[each_atom][2]) + '\n'
                        if atomic_numbers[each_atom] == 6:
                            init_xyz_string = init_xyz_string + 'C' + ' ' + str(xyz_positions[each_atom][0]) + ' ' + str(xyz_positions[each_atom][1]) + ' ' + str(xyz_positions[each_atom][2]) + '\n'
                        if atomic_numbers[each_atom] == 7:
                            init_xyz_string = init_xyz_string + 'N' + ' ' + str(xyz_positions[each_atom][0]) + ' ' + str(xyz_positions[each_atom][1]) + ' ' + str(xyz_positions[each_atom][2]) + '\n'
                        if atomic_numbers[each_atom] == 8:
                            init_xyz_string = init_xyz_string + 'O' + ' ' + str(xyz_positions[each_atom][0]) + ' ' + str(xyz_positions[each_atom][1]) + ' ' + str(xyz_positions[each_atom][2]) + '\n'
                        if atomic_numbers[each_atom] == 9:
                            init_xyz_string = init_xyz_string + 'F' + ' ' + str(xyz_positions[each_atom][0]) + ' ' + str(xyz_positions[each_atom][1]) + ' ' + str(xyz_positions[each_atom][2]) + '\n'

#                    scratch_xyz_filename = self.scratch_file + '.xyz' 
                    self.init_dataset_xyzfile.write(init_xyz_string)

                    #make transformed xyz
                    trans_xyz_string = ''
                    countH_inxyz = 0

                    for this_line_index, this_line in enumerate(init_xyz_string.split('\n')):
                        if this_line_index == 0:
                            this_line = this_line.replace(this_line,'\n'+str(int(this_line)-2)+'\n')

                        if 'H' in this_line:
                            if countH_inxyz == H1_index_in_xyzmol[0]:
                                this_line = this_line.replace(this_line,'')
                            if countH_inxyz == H2_index_in_xyzmol[0]:
                                this_line = this_line.replace(this_line,'')
                            countH_inxyz = countH_inxyz + 1
                    
                        trans_xyz_string = trans_xyz_string + this_line +'\n'
                    
                    print('before',trans_xyz_string)

                    #remove empty lines
                    lines = trans_xyz_string.split('\n')
                    trans_xyz_nonempty_lines = [line for line in lines if line.strip()]

                    print('nonempty',trans_xyz_nonempty_lines)

                    # Join the non-empty lines back into a string
                    trans_xyz_string2 = '\n'.join(trans_xyz_nonempty_lines)+'\n'

                    print('after',trans_xyz_string2)
                    self.trans_dataset_xyzfile.write(trans_xyz_string2)

                else:
                    pass

        self.init_dataset_xyzfile.close()
        self.trans_dataset_xyzfile.close()
        
        print('init_xyz_string',init_xyz_string)
        print('trans_xyz_string',trans_xyz_string)


        utils.generate(self.init_dataset_xyzfilepath,self.available_properties,self.n_molecules[1])
        utils.generate(self.trans_dataset_xyzfilepath,self.available_properties,self.n_molecules[1])
