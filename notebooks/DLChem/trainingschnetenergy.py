# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 23:09:30 2021

@author: aelsamma
"""

import os

import torch

import schnetpack.nn 
import schnetpack.data
import schnetpack as spk
from schnetpack.datasets import QM9


import schnetpack.train as trn
from torch.optim import Adam

import numpy as np
import matplotlib.pyplot as plt
from ase.units import kcal, mol


def mse_loss(batch, result):
    diff = batch[QM9.U0]-result[QM9.U0]
    err_sq = torch.mean(diff ** 2)
    return err_sq
    
def train_schnet_model(trained_model_path,qm9data,num_train,num_val,n_in,n_atomic_basis,n_filters,n_gaussians,n_interactions,cutoff,device):
    trained_model_path = trained_model_path
    if not os.path.exists(trained_model_path):
        os.makedirs(trained_model_path)

    qm9data=qm9data
    train, val, test = spk.train_test_split(
        data=qm9data,
        num_train=num_train,
        num_val=num_val,
        split_file=os.path.join(trained_model_path,'split.npz')
        )


    train_loader = spk.AtomsLoader(train, batch_size=100, shuffle=True)
    val_loader = spk.AtomsLoader(val, batch_size=100)


    atomrefs = qm9data.get_atomref(QM9.U0)
    print('U0 of hyrogen:', '{:.2f}'.format(atomrefs[QM9.U0][1][0]), 'eV')
    print('U0 of carbon:', '{:.2f}'.format(atomrefs[QM9.U0][6][0]), 'eV')
    print('U0 of oxygen:', '{:.2f}'.format(atomrefs[QM9.U0][8][0]), 'eV')
    print('U0 of nitrogen:', '{:.2f}'.format(atomrefs[QM9.U0][7][0]), 'eV')

    means, stddevs = train_loader.get_statistics(
        QM9.U0, divide_by_atoms=True, single_atom_ref=atomrefs)


    print('Mean atomization energy / atom:', means[QM9.U0])
    print('Std. dev. atomization energy / atom:', stddevs[QM9.U0])

    # Define output model and property to be predicted
    schnet = spk.representation.SchNet(
        n_atom_basis=30, n_filters=30, n_gaussians=20, n_interactions=3,
        cutoff=4., cutoff_network=spk.nn.cutoff.CosineCutoff
    )
    output_U0 = spk.atomistic.Atomwise(n_in=n_in, atomref=atomrefs[QM9.U0], property=QM9.U0,
                                   mean=means[QM9.U0], contributions=True, stddev=stddevs[QM9.U0])
    model = spk.AtomisticModel(representation=schnet, output_modules=output_U0)

    # loss function


    # build optimizer
    optimizer = Adam(model.parameters(), lr=1e-2)

    loss = trn.build_mse_loss([QM9.U0])

    metrics = [spk.metrics.MeanAbsoluteError(QM9.U0)]
    hooks = [
        trn.CSVHook(log_path=trained_model_path, metrics=metrics),
        trn.ReduceLROnPlateauHook(
            optimizer,
            patience=5, factor=0.8, min_lr=1e-6,
            stop_after_min=True
        )
    ]

    trainer = trn.Trainer(
        model_path=trained_model_path,
        model=model,
        hooks=hooks,
        loss_fn=loss,
        optimizer=optimizer,
        keep_n_checkpoints = 10,
        train_loader=train_loader,
        validation_loader=val_loader,
    )

    device = device # change to 'cpu' if gpu is not available
    n_epochs = 1000 # takes about 10 min on a notebook GPU. reduces for playing around
    trainer.train(device=device, n_epochs=n_epochs)


    results = np.loadtxt(os.path.join(trained_model_path, 'log.csv'), skiprows=1, delimiter=',')

    time = results[:,0]-results[0,0]
    learning_rate = results[:,1]
    train_loss = results[:,2]
    val_loss = results[:,3]
    val_mae = results[:,4]

    print('Final validation MAE:', np.round(val_mae[-1], 2), 'eV =',
          np.round(val_mae[-1] / (kcal/mol), 2), 'kcal/mol')

    plt.figure(figsize=(14,5))
    plt.subplot(1,2,1)
    plt.plot(time, val_loss, label='Validation')
    plt.plot(time, train_loss, label='Train')
    plt.yscale('log')
    plt.ylabel('Loss [eV]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.subplot(1,2,2)
    plt.plot(time, val_mae)
    plt.ylabel('mean abs. error [eV]')
    plt.xlabel('Time [s]')
    plt.show()

    checkpoint_path = trained_model_path + '/trained.pth'
    torch.save(model.state_dict(), checkpoint_path)
    

def use_model(idx,model_file_dir,qm9data):
     ### #Loads the best_model after schnet training.
    #user evaluates the model on the QM9 database
    model_file = model_file_dir + "best_model"
    split_file = model_file_dir + 'split.npz'
    idx=idx
    # Load best model where it was stored after training
    
    model = torch.load(model_file, map_location=torch.device('cpu'))

    # Load split file 
    train, val, test = spk.data.train_test_split(qm9data, split_file=split_file)
    
    at, props = qm9data.get_properties(idx)
    
    device = 'cpu'
    
    # load data for molecule
    
    # load spk calculator
    calculator = spk.interfaces.SpkCalculator(model=model, device=device, energy=QM9.mu)
    # set calculator on molecule
    at.set_calculator(calculator)
    
    converter = spk.data.AtomsConverter(device=device)
    # convert qm9 data to machine-readable form
    inputs = converter(at)
    
    positions = props['_positions']
    positions = positions.detach().numpy()
    
    for i in range(len(props['_atomic_numbers'])):
        if props['_atomic_numbers'][i] == 1:
            print('H ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]))
        if props['_atomic_numbers'][i] == 6:
            print('C ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]))
        if props['_atomic_numbers'][i] == 7:
            print('N ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]))
        if props['_atomic_numbers'][i] == 8:
            print('O ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]))
        if props['_atomic_numbers'][i] == 9:
            print('F ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]))
                        
    pred = model(inputs)
    
    print('Prediction:', pred[QM9.U0].detach().cpu().numpy()[0,0])
    print('Keys:', list(inputs.keys()))
    print('Truth:', props[QM9.U0].cpu().numpy()[0])