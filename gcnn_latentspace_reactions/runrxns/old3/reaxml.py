from schnetpack.datasets import QM9
from schnetpack import AtomsData

from numpy import genfromtxt
import numpy as np

from numpy import savetxt

class rxns():
    def __init__(self, db_filepath, qm9, available_properties,mol_range,label_filepath,label1,timesteps,scale_factor,output_filepath,initial_H_atom_idx) -> None:

        #initialize molecular database parameters
        self.db_filepath = db_filepath
        self.qm9 = qm9
        self.mol_range = mol_range
        self.available_properties = available_properties
        self.label_filepath = label_filepath
        self.label1 = label1
        self.timesteps = timesteps
        self.scale_factor = scale_factor
        self.output_filepath = output_filepath
        self.initial_H_atom_idx = initial_H_atom_idx

        #load db_data (molecular database file)
        if qm9 == True:
            db_data = QM9(db_filepath,download=False,remove_uncharacterized=True)
        else:
            db_data = AtomsData(db_filepath,available_properties=self.available_properties)
    
        self.db_data = db_data

        #load lda label column from the label file
        self.labeldata = genfromtxt(self.label_filepath,delimiter=' ',encoding='utf-8-sig')

        self.output_file = open(self.output_filepath,mode='w')

        self.atomic_number_dict = {1: "H",6: "C",7: "N",8: "O",9: "F"}

    def oxidize(self):

        #important to initialize a total atom count to keep track of the
        #number of atoms in each molecule as we extract labels from the embs file
        total_H_atom_count = self.initial_H_atom_idx
            
        #count molecules to keep track of how many satisfy fg label for reaction
        molecule_count = 0
        for each_molecule in range(self.mol_range[0],self.mol_range[1]):



            at, props = self.db_data.get_properties(each_molecule)


            xyz_coords = props['_positions'].detach().numpy()
            atomic_numbers = props['_atomic_numbers'].detach().numpy()
            number_atoms = len(atomic_numbers)

            number_O = np.count_nonzero(atomic_numbers==8)
            if number_O == 1:
                for each_atom in range(number_atoms):
                    
                    #perform a shifting operation only for Hydrogen-type OH-CH2 type atoms
                    #note that you will have to only choose one of the two CH2 hydrogens
                    #to shift for the oxidation reaction
                    if props['_atomic_numbers'][each_atom] == 1:

                        #get the label from the H-atoms autolabel file
                        label = self.labeldata[total_H_atom_count][0]


                        if label == self.label1:
                            for timestep in range(self.timesteps):
                                
                                molecule_count = molecule_count + 1
                                #find bond vector, so you can dissociate H along the bond direction/axis
                                #the easiest way to do this is to find that oxygen atom that is closest
                                #take the two vectors and calculate a vector difference
                                O_nbr_idx_dists = []
                                H_nbr_idx_dists = []
                                for each_neighbor in range(number_atoms):
                                    if atomic_numbers[each_neighbor] == 8:
                                        O_nbr_distance = np.linalg.norm(xyz_coords[each_atom]-xyz_coords[each_neighbor])
                                        O_nbr_idx_dists.append([each_neighbor,O_nbr_distance])
                                    if atomic_numbers[each_neighbor] == 1:
                                        H_nbr_distance = np.linalg.norm(xyz_coords[each_atom]-xyz_coords[each_neighbor])
                                        if H_nbr_distance > 0:
                                            H_nbr_idx_dists.append([each_neighbor,H_nbr_distance])

                                #find the closest oxygen neighbor
                                O_nbr_idx_dists = np.array(O_nbr_idx_dists)
                                O_min_idx = np.argmin(O_nbr_idx_dists[:,1])
        
                                #find the closest hydrogen neighbor
                                H_nbr_idx_dists = np.array(H_nbr_idx_dists)
                                H_min_idx = H_nbr_idx_dists[np.argmin(H_nbr_idx_dists[:,1])][0]
                                
                                #calculate vector diff (bond axis)
                                if timestep == 0:
                                    bond_vector =   xyz_coords[each_atom] - xyz_coords[O_min_idx]
                                else:
                                    bond_vector = bond_vector

                                new_xyz_coords = [xyz_coords[each_atom][dim] + timestep*self.scale_factor*bond_vector[dim] for dim in range(3)]
                                nearby_H_new_xyz_coords = [xyz_coords[int(H_min_idx)][dim] + timestep*self.scale_factor*bond_vector[dim] for dim in range(3)]

                                xyz_coords[each_atom] = new_xyz_coords
                                xyz_coords[int(H_min_idx)] = nearby_H_new_xyz_coords

                                self.output_file.write(str(number_atoms)+'\n' + '0.00000' + '\n')
                                for loop_atom in range(number_atoms):                   
                                    element_of_atom = self.atomic_number_dict.get(atomic_numbers[loop_atom])
                                    self.output_file.write(element_of_atom + ' ' +  str(xyz_coords[loop_atom][0]) + ' ' +  str(xyz_coords[loop_atom][1]) + ' ' +  str(xyz_coords[loop_atom][2]) + '\n')
                
                        total_H_atom_count = total_H_atom_count + 1

            else:               
                #increase H atom count 
                total_H_atom_count = total_H_atom_count + np.count_nonzero(atomic_numbers==1)
            
        self.output_file.close()
      
#        with open(self.output_filepath.replace('.txt','mol_count.txt'), 'w') as file:
#            file.write(str(molecule_count))
            
        print(molecule_count)


        

