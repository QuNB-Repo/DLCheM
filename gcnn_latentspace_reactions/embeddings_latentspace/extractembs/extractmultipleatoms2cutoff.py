from schnetpack.datasets import QM9
import schnetpack as spk
import torch
import os
from schnetpack import AtomsData

import numpy as np
from numpy import savetxt, genfromtxt

def extract(db_file_path,model_filepath,save_filepath,layer,n_features,n_molecules,n_atoms,cutoff,target_element):

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

        number_atoms = len(positions)
        

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