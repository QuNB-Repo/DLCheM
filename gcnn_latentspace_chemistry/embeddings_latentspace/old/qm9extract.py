'''
This is the main tool. I have tried defining it into
its own function, but loading model checkpoints does 
not seem to allow that. Nonetheless, it can be easily
controlled from here

This extracts the inner layers of a trained model, 
first by loading saved checkpoint of model
then running a forward pass through the model

As soon as this is placed inside a function,
the ability to load the model is lost for some reason.

'''

import torch
import schnetpack as spk
import math
from schnetpack.datasets import QM9
import os

import numpy as np   
from numpy import savetxt

import pandas as pd

from qm9input import *

#read in label file in order to output specific functional groups
#label = pd.read_csv(label_file,delimiter=',')

qm9data = QM9(qm9_file, download=False, remove_uncharacterized=True)

def hook(self, inp_tensor, out_tensor):
    # Self is included and refers to the model class
    # Global allows us to utilize embedding_output outside the current function scope
    global layer
    #Update the embedding_output variable to be equal to our output tensor
    layer=out_tensor 

# Load split file 
train, val, test = spk.data.train_test_split(qm9data,split_file=split_file)

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
load_checkpoint = torch.load(checkpoint_path,map_location=torch.device('cpu'))


#qm9_i6_30f_20g-1000-500-4_300.pth
# load model's state dictionary from saved checkpoint
model.load_state_dict(load_checkpoint)


#set up device for forward pass
device='cpu'

# load atoms converter 
converter = spk.data.AtomsConverter(device=device)

dataH = np.zeros((1,n_filters))
dataaeH = np.zeros((1,1))
dataC = np.zeros((1,n_filters))
dataaeC = np.zeros((1,1))
dataN = np.zeros((1,n_filters))
dataaeN = np.zeros((1,1))
dataO = np.zeros((1,n_filters))
dataaeO = np.zeros((1,1))
dataF = np.zeros((1,n_filters))
dataaeF = np.zeros((1,1))
dataall = np.zeros((1,n_filters))
dataallae = np.zeros((1,1))
countO = 0
label  = np.zeros((1,1))
for idx in range(start,end):
    print(idx)

    at, props = qm9data.get_properties(idx)
    inputs = converter(at)
    
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
    
    
    rep = emb + int0 +int1 +int2
    
    number_atoms = len(props['_positions'])
    

    from schnetpack.atomistic.output_modules import yi
    

    print('extract rep', rep)

    rows = np.zeros((number_atoms,n_filters))
    for i in range(number_atoms):
        for j in range(n_filters):
            rows[i][j] = rep[0][i][j]

#    print(rows)
    yi=yi.detach().numpy()
    #save the vector of every oxygen atom encountered
    
    from schnetpack.representation.schnet import xemb

    print('schnet x', xemb)
    
#    print(idx)
    for i in range(number_atoms):
#        if props['_atomic_numbers'][i] == 8:
#            if label['Target'][countO] == 2:
#                print(idx)
#                datao = np.vstack((datao,rows[i]))
#                dataoae = np.vstack((dataoae,yi[0][i])) 
#            countO = countO + 1

        dataall = np.vstack((dataall,rows[i])) 
        dataallae = np.vstack((dataallae,yi[0][i]))
        
        if props['_atomic_numbers'][i] == 1:
            print('H ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            dataH = np.vstack((dataH,rows[i])) 
            dataaeH = np.vstack((dataaeH,yi[0][i]))
            label = np.vstack((label,0)) 
        elif props['_atomic_numbers'][i] == 8:
            print('O ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            dataO = np.vstack((dataO,rows[i])) 
            dataaeO = np.vstack((dataaeO,yi[0][i]))
            label = np.vstack((label,3)) 
        elif props['_atomic_numbers'][i] == 7:
            print('N ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            dataN = np.vstack((dataN,rows[i]))    
            dataaeN = np.vstack((dataaeN,yi[0][i]))
            label = np.vstack((label,2)) 
        elif props['_atomic_numbers'][i] == 6:
            print('C ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            dataC = np.vstack((dataC,rows[i]))    
            dataaeC = np.vstack((dataaeC,yi[0][i]))
            label = np.vstack((label,1)) 
        elif props['_atomic_numbers'][i] == 9:
            print('F ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            dataF = np.vstack((dataF,rows[i]))    
            dataaeF = np.vstack((dataaeF,yi[0][i]))
            label = np.vstack((label,4)) 
    
dataH = np.delete(dataH, 0 ,axis=0)
dataaeH = np.delete(dataaeH, 0 ,axis=0)

dataC = np.delete(dataC, 0 ,axis=0)
dataaeC = np.delete(dataaeC, 0 ,axis=0)

dataN = np.delete(dataN, 0 ,axis=0)
dataaeN = np.delete(dataaeN, 0 ,axis=0)

dataO = np.delete(dataO, 0 ,axis=0)
dataaeO = np.delete(dataaeO, 0 ,axis=0)

dataF = np.delete(dataF, 0 ,axis=0)
dataaeF = np.delete(dataaeF, 0 ,axis=0)

dataall = np.delete(dataall, 0 ,axis=0)
dataallae = np.delete(dataallae, 0 ,axis=0)

savetxt(save_filepathall, dataall, delimiter=',',encoding='utf-8-sig')
save_filepathallae = save_filepathall.replace('.csv','ae.csv')
savetxt(save_filepathallae,dataallae,delimiter=',',encoding='utf-8-sig')

savetxt(save_filepathH,dataH,delimiter=',',encoding='utf-8-sig')
save_filepathHae = save_filepathH.replace('.csv','ae.csv')
savetxt(save_filepathHae,dataaeH,delimiter=',',encoding='utf-8-sig')

savetxt(save_filepathC,dataC,delimiter=',',encoding='utf-8-sig')
save_filepathCae = save_filepathC.replace('.csv','ae.csv')
savetxt(save_filepathCae,dataaeC,delimiter=',',encoding='utf-8-sig')

savetxt(save_filepathN,dataN,delimiter=',',encoding='utf-8-sig')
save_filepathNae = save_filepathN.replace('.csv','ae.csv')
savetxt(save_filepathNae,dataaeN,delimiter=',',encoding='utf-8-sig')

savetxt(save_filepathO,dataO,delimiter=',',encoding='utf-8-sig')
save_filepathOae = save_filepathO.replace('.csv','ae.csv')
savetxt(save_filepathOae,dataaeO,delimiter=',',encoding='utf-8-sig')

savetxt(save_filepathF,dataF,delimiter=',',encoding='utf-8-sig')
save_filepathFae = save_filepathF.replace('.csv','ae.csv')
savetxt(save_filepathFae,dataaeF,delimiter=',',encoding='utf-8-sig')

savetxt(save_all_label,label,delimiter=',',encoding='utf-8-sig')
