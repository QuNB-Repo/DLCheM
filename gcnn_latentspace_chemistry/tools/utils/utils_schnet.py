from schnetpack.datasets import QM9

import torch

import schnetpack.nn 
import schnetpack.data
import schnetpack as spk

from schnetpack import AtomsData

import schnetpack.train as trn
from torch.optim import Adam

import numpy as np
import matplotlib.pyplot as plt
from ase.units import kcal, mol

import os


def run_model(dataset_filepath,model_filepath,idx,qm9):
#Load testing data

    if qm9 == True:
        qm9data = QM9(dataset_filepath,download=False,
                    remove_uncharacterized=True)
    else:
        qm9data = AtomsData(dataset_filepath,available_properties=['energy'])

#    print(dataset_filepath)
#    print(model_filepath)

    model = torch.load(model_filepath, map_location=torch.device('cpu'))

    at, props = qm9data.get_properties(idx)

    positions = props['_positions']

    #set device
    device = 'cpu'

    # load the data converter which will convert the data to machine-readable format for the algorithm
    converter = spk.data.AtomsConverter(device=device)
    # convert data to machine-readable form
    inputs = converter(at)

    #use model on inputs  
    pred = model(inputs)

#    print('Keys:', list(inputs.keys()))
#    print('Prediction:', pred[QM9.U0].detach().cpu().numpy()[0,0])
#    print('Truth:', props[QM9.U0][0])
    return pred[QM9.U0].detach().cpu().numpy()[0,0]

#this below is important for training (maybe make a class next time when you learn how to train other properties, the class controls how to train each, and defines MSE loss separately to be used by each)
def mse_loss(batch, result):
    diff = batch[QM9.U0]-result[QM9.U0]
    err_sq = torch.mean(diff ** 2)
    return err_sq


def train_energy(trained_model_path,qm9,dataset_filepath,num_train,num_val,n_in,n_atom_basis,n_filters,n_gaussians,n_interactions,cutoff,lr,batch_size,device):
    trained_model_path = trained_model_path
    if not os.path.exists(trained_model_path):
        os.makedirs(trained_model_path)

    if qm9 == True:
        qm9data = QM9(dataset_filepath,download=False,
                    remove_uncharacterized=True)
    else:
        qm9data = AtomsData(dataset_filepath,available_properties=['energy'])

    train, val, test = spk.train_test_split(
        data=qm9data,
        num_train=num_train,
        num_val=num_val,
        split_file=os.path.join(trained_model_path,'split.npz')
        )


    train_loader = spk.AtomsLoader(train, batch_size=batch_size, shuffle=True)
    val_loader = spk.AtomsLoader(val, batch_size=batch_size)


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
        n_atom_basis=n_atom_basis, n_filters=n_filters, n_gaussians=n_gaussians, n_interactions=n_interactions,
        cutoff=cutoff, cutoff_network=spk.nn.cutoff.CosineCutoff
    )
    output_U0 = spk.atomistic.Atomwise(n_in=n_in, atomref=atomrefs[QM9.U0], property=QM9.U0,
                                   mean=means[QM9.U0], contributions=True, stddev=stddevs[QM9.U0])
    model = spk.AtomisticModel(representation=schnet, output_modules=output_U0)

    # loss function


    # build optimizer
    optimizer = Adam(model.parameters(), lr=lr)

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
    

def extra():
    print(len(positions))

    from schnetpack.atomistic.output_modules import yi

    for i in range(len(positions)):
        if props['_atomic_numbers'][i] == 1:
            print('H ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            print(yi[0][i])
        if props['_atomic_numbers'][i] == 6:
            print('C ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            print(yi[0][i])
        if props['_atomic_numbers'][i] == 7:
            print('N ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            print(yi[0][i])
        if props['_atomic_numbers'][i] == 8:
            print('O ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            print(yi[0][i])
        if props['_atomic_numbers'][i] == 9:
            print('F ' + str(props['_positions'][i][0]) + ' ' + str(props['_positions'][i][1]) + ' ' + str(props['_positions'][i][2]))
            print(yi[0][i])
