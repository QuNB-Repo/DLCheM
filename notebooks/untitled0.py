#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 07:21:57 2021

@author: amerelsamman
"""


import schnetpack as spk
import torch
import schnetpack.nn 
import schnetpack.data
import pandas as pd
import scipy.linalg as la
from schnetpack.datasets import QM9

#Intermediare embedding output pipeline

#Load last checkpoint file from from trained model
Checkpoint = './trained_models/qm911/checkpoints/checkpoint-5.pth.tar'

#Load the Best model where it is stored after training
best_model = torch.load('./trained_models/qm911/best_model')
model = best_model

model.load_state_dict(torch.load(Checkpoint),strict=False)

# Download QM9 dataset to use in evaluating model
qm9data = QM9('./qm9.db', download=True, remove_uncharacterized=True)

# Load split file for the train, validation and test 
train, val, test = spk.data.train_test_split(qm9data, split_file='./trained_models/qm911/split.npz')



#set up device and atoms converter for input
device = 'cpu'

converter = spk.data.AtomsConverter(device=device)

test_loader = spk.AtomsLoader(test, batch_size=100)
converter = spk.data.AtomsConverter(device=device)

at, props = qm9data.get_properties(idx=1)

calculator = spk.interfaces.SpkCalculator(model=best_model, device=device, energy=QM9.U0)
at.set_calculator(calculator)


inputs = converter(at)

print('Prediction:', at.get_total_energy())


#We choose None to instatiate the variable originally
embedding_output=None
#Define a hooking function that will fetch the out_tensor 
#from the layer of interest given an input tensor
def embedding_hook(self, inp_tensor, out_tensor):
    # Self is included and refers to the model class
    # Global allows us to utilize embedding_output outside the current function scope
    global embedding_output
    #Update the embedding_output variable to be equal to our output tensor
    embedding_output=out_tensor

#Hook our function to the embedding layer during the forward pass
model.representation.embedding.register_forward_hook(embedding_hook)

#Check if None
print(embedding_output)

#Forward pass of the tensor inputs
model(inputs)

#embedding_output SHOULD be changed
print(embedding_output)




interactions_output=None

def interaction_hook(self, inp_tensor, out_tensor):
    global interactions_output
    interactions_output=out_tensor

model.representation.interactions.register_forward_hook(interaction_hook)

#Check if None
print(interactions_output)

#Forward pass of the tensor inputs
model(inputs)

#embedding_output SHOULD be changed
print(interactions_output)
