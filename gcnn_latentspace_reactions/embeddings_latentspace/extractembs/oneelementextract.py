import torch
import schnetpack as spk
import math
from schnetpack import AtomsData
from schnetpack.datasets import QM9
import os

import numpy as np
from numpy import savetxt, genfromtxt

def extract(db_file_path,model_filepath,save_filepath,start,end,atom_num,label_id,use_label,label_filepath,newlabel,layer,n_features):
    #dataset used for testing embedding vectors
    dataset = AtomsData(db_file_path, available_properties=['energy'])

    if use_label == True:
        label = genfromtxt(label_filepath,delimiter=' ',encoding='utf-8-sig')

    model = torch.load(model_filepath, map_location=torch.device('cpu'))

    #set device
    device = 'cpu'

    # load the data converter which will convert the data to machine-readable format for the algorithm
    converter = spk.data.AtomsConverter(device=device)

    data = open(save_filepath,mode='a',encoding='utf-8-sig')
    save_filepathae = save_filepath.replace('.csv','ae.csv')
    dataae = open(save_filepathae,mode='a',encoding='utf-8-sig')

    count_atom = 0
    for idx in range(start,end):
        print(idx)

        #get molecule properties
        at, props = dataset.get_properties(idx)

        #convert to machine readable format
        inputs = converter(at)

        #use model on inputs  
        pred = model(inputs)

        #get atomic goordinates
        positions = props['_positions']

        from schnetpack.representation.schnet import xs
        from schnetpack.atomistic.output_modules import yi

        x =  xs[layer].cpu().detach().numpy()
        yi  = yi.detach().numpy()
        
        n_atoms = len(positions)

        for i in range(n_atoms):
            if props['_atomic_numbers'][i] == atom_num:
                if use_label == True:
                    if label[count_atom][0] == label_id:
                        row = ''
                        for k in range(n_features):
                            row = row + str(x[0][i][k]) + ', '
                        row = row + '\n'

                        rowae = '' + str(yi[0][i][0]) + '\n'
                        
                        data.write(row)
                        dataae.write(rowae)

                else:
                    row = ''
                    for k in range(n_features):
                        row = row + str(x[0][i][k]) + ', '
                    row = row + '\n'

                    rowae = '' + str(yi[0][i][0]) + '\n'
                    
                    data.write(row)
                    dataae.write(rowae)
                count_atom = count_atom + 1

    if use_label == True:
        new_label = np.delete(new_label,0,axis=0)
        new_label_filepath = save_filepath.replace('.csv','label.csv')
        savetxt(new_label_filepath,new_label,delimiter=',') 
    
    data.close()
    dataae.close()