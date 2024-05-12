import torch
import schnetpack as spk
import math
from schnetpack import AtomsData
from schnetpack.datasets import QM9
import os

import numpy as np
from numpy import savetxt, genfromtxt

from input import *

#qm9 data required to load checkpoint for extraction
qm9data = QM9(qm9_file, download=False, remove_uncharacterized=True)

#dataset used for testing embedding vectors
dataset = AtomsData(db_file_path, available_properties=['energy'])


label = genfromtxt(label_filepath,delimiter=',',encoding='utf-8-sig')


# define schnet variables as they are saved in checkpoint (to load checkpoint and extract internal layers, this is required)

def hook(self, inp_tensor, out_tensor):
    # Self is included and refers to the model class
    # Global allows us to utilize embedding_output outside the current function scope
    global layer
    #Update the embedding_output variable to be equal to our output tensor
    layer=out_tensor 

vecs = np.zeros((1,31))
count_atom = 0

for idx in range(index):
    # Load split file 
#    train, val, test = spk.data.train_test_split(qm9data,split_file=split_file)

    # Load atom ref data 
    atomrefs = qm9data.get_atomref(QM9.U0)

    # Define SchNet representation model

    schnet = spk.representation.SchNet(
    n_atom_basis=n_atom_basis, n_filters=n_filters, n_gaussians=n_gaussians, n_interactions=n_interactions,
    cutoff=cutoff , cutoff_network=spk.nn.cutoff.CosineCutoff
    )

    # Define SchNet output model and property to be predicted

    output_U0 = spk.atomistic.Atomwise(n_in=n_filters, atomref=atomrefs[QM9.U0])

    # Define atomistic model

    model = spk.AtomisticModel(representation=schnet,output_modules=output_U0)

    # Load saved checkpoint file
    load_checkpoint = torch.load(checkpoint_path, map_location=torch.device('cpu'))


    #qm9_i6_30f_20g-1000-500-4_300.pth
    # load model's state dictionary from saved checkpoint
    model.load_state_dict(load_checkpoint)


    #set up device for forward pass
    device='cpu'

    # load spk calculator
    calculator = spk.interfaces.SpkCalculator(model=model, device=device, energy=QM9.U0)
    converter = spk.data.AtomsConverter(device=device)

    at, props = dataset.get_properties(idx)
    inputs = converter(at)
    number_atoms = len(props['_atomic_numbers'])

    print(at)

    layer = None

    model.representation.embedding.register_forward_hook(hook)
    model(inputs)

    emb = layer.clone()
    emb = layer.detach().numpy()

    layer = None

    model.representation.interactions[0].register_forward_hook(hook)

    model(inputs)
    int0 = layer.clone()
    int0 = int0.detach().numpy()    

    layer = None

    model.representation.interactions[1].register_forward_hook(hook)

    model(inputs)
    int1 = layer.clone()
    int1 = int1.detach().numpy()   

    layer = None 

    model.representation.interactions[2].register_forward_hook(hook)

    model(inputs)   

    int2 = layer.clone()
    int2 = int2.detach().numpy()   

    rep = emb+int0+int1+int2

#    from schnetpack.atomistic.output_modules import yi

#    yi=yi.detach().numpy()
    
    rows = np.zeros((number_atoms,31))
    for i in range(number_atoms):
        for j in range(30):
            rows[i][j] = rep[0][i][j]
#        rows[i][30] = yi[0][i]

    #should have a counter for where you are in the labelfile! 

    for i in range(number_atoms):
        if props['_atomic_numbers'][i] == atom_num:
            if use_label == True:
                if label[count_atom][0] == label_id:
                    vecs = np.vstack((vecs,rows[i]))

            else:
                vecs = np.vstack((vecs,rows[i]))
            count_atom = count_atom + 1
#        if props['_atomic_numbers'][i] == 7:
#            vecs = np.vstack((vecs,rows[i]))
#        if props['_atomic_numbers'][i] == 6 and 8 not in props['_atomic_numbers'] and 7 not in props['_atomic_numbers']: 
#            vecs = np.vstack((vecs,rows[i]))

vecs = np.delete(vecs, 0, axis=0)
savetxt(save_filepath,vecs,delimiter=',')

lab = np.delete(lab,0,axis=0)
lab_filepath = save_filepath.replace('.csv','label.csv')
savetxt(lab_filepath,lab,delimiter=',') 