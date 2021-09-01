# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 18:13:58 2021

@author: aelsamma
"""

import torch
from schnetpack.datasets import QM9
import schnetpack.nn 
import schnetpack.data
import schnetpack as spk


# define the train/val split file which is in the folder where your trained model was saved
split_file = '../../../data/trainedmodels/qm9energy10000-30/split.npz'
# define your best model file, this defines what your model has learned and it is found in your trained folder
model_file = '../../../data/trainedmodels/qm9energy10000-30/best_model'
#load it with torch, and use CPU if you trained with CPU
model = torch.load(model_file, map_location=torch.device('cpu'))

#load your data (I used QM9, but you can load any with AtomsData, has to be db file)
qm9data = QM9('../../../data/datasets/QM9/qm9.db')

# Load split file 
train, val, test = spk.data.train_test_split(qm9data, split_file=split_file)

#choose a datapoint from your data set to test and get its "at" and "props"
idx=0
at, props = qm9data.get_properties(idx)

#set device
device = 'cpu'

# load the data converter which will convert the data to machine-readable format for the algorithm
converter = spk.data.AtomsConverter(device=device)
# convert data to machine-readable form
inputs = converter(at)

positions = props['_positions']
positions = positions.detach().numpy()

#use model on inputs  
pred = model(inputs)

print('Keys:', list(inputs.keys()))
print('Prediction:', pred[QM9.U0].detach().cpu().numpy()[0,0])
print('Truth:', props[QM9.U0].cpu().numpy()[0])