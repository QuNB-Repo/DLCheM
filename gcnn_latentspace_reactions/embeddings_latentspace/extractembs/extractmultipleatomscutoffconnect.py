from schnetpack.datasets import QM9
import schnetpack as spk
import torch
import os
from schnetpack import AtomsData
from extractembeddings.utils import utils_connect

import numpy as np
from numpy import savetxt, genfromtxt

def extract(db_file_path,model_filepath,save_filepath,layer,n_features,n_molecules,n_atoms,cutoff,target_element,label_dir):
    

    dataset = AtomsData(db_file_path, available_properties=['energy'])

    device = 'cpu'

    model = torch.load(model_filepath, map_location=torch.device(device))

    converter = spk.data.AtomsConverter(device=device)



    data = open(save_filepath,mode='a',encoding='utf-8-sig')
    save_filepathae = save_filepath.replace('.csv','ae.csv')
    dataae = open(save_filepathae,mode='a',encoding='utf-8-sig')
    
    for idx in range(120,n_molecules):

        if idx % 1000 == 0:
            print(idx)
        
        at, props = dataset.get_properties(idx)


        inputs = converter(at)

        model(inputs)
        
        positions = props['_positions']
        num_atoms = len(props['_atomic_numbers'])


        #find index of target element
        for i in range(num_atoms):
            if props['_atomic_numbers'][i] == target_element:
                target_index = i

        num_carbons = 0
        num_hydrogens = 0
        dist_tooxygen = []
        for i in range(num_atoms):
            if props['_atomic_numbers'][i] == 6:
                num_carbons = num_carbons + 1
            if props['_atomic_numbers'][i] == 1:
                num_hydrogens = num_hydrogens + 1

            #calculate distance to target index and add to dist_tooxygen...
            distance = np.linalg.norm(props['_positions'][i,:]-props['_positions'][target_index,:])

            dist_tooxygen.append(distance)

        from schnetpack.representation.schnet import xs
        from schnetpack.atomistic.output_modules import yi

        x = xs[layer].cpu().detach().numpy()
        yi  = yi.detach().numpy()
        #argsort embedding vector data according to dist_tooxygen vector
        dist_tooxygen = np.transpose(np.array(dist_tooxygen))

#        NOTE DO NOT ARGSORT  ... WHAT if it needs to be ordered according to connectivity      
        x[0,:] = x[0][dist_tooxygen.argsort()]
        props['_atomic_numbers'] = props['_atomic_numbers'][dist_tooxygen.argsort()]

        number_atoms = len(props['_atomic_numbers'])
        
        #sort the embedding vectors according to connectivity of each starting with target atom
        #write xyz from db
        name_xyz = utils_connect.write_xyz_from_db(props,idx,label_dir)
        #convert xyz to mol
        name_mol = utils_connect.xyz_to_mol(idx,label_dir)

        #construct total_branches_list
        connect_list, element_list = utils_connect.connections_list(name_mol,number_atoms)

        print(connect_list)
        print(element_list)
        #order connection matrix using some standard
        n_branches=len(connect_list)
        totalconnect_list, totalelement_list = utils_connect.order_branches(connect_list, element_list,n_branches)

        print('ordconnect',totalconnect_list)
        print('ordelemen',totalelement_list)

        break
        if number_atoms >= n_atoms:
            countatoms = 0 
            row = '' 
            for i in range(number_atoms):
                if countatoms < cutoff:
                    countatoms = countatoms + 1

                    for k in range(n_features):
                        row = row + str(x[0][i][k]) + ', '
            
                rowae = str(yi[0][i]) 

            row = row + ' \n'
            rowae = rowae + ' \n'
            data.write(row)
            dataae.write(rowae)

    data.close()
    dataae.close()