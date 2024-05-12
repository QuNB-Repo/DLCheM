from schnetpack.datasets import QM9
import schnetpack as spk
import torch 

import numpy as np
from numpy import genfromtxt, savetxt

from schnetpack import AtomsData


from oldtools.utils import utils_xyzdbmol

def extract(data_filepath,qm9,available_properties,model_filepath,n_features,start,end,save_filepath,layer,labelofilepath,pkafilepath,labelnfilepath=''):
    #Load testing data
    if qm9 == True:
        dataset = QM9(data_filepath,download=False,
                    remove_uncharacterized=True)
    else:
        dataset = AtomsData(data_filepath,available_properties=available_properties)

    model = torch.load(model_filepath, map_location=torch.device('cpu'))

    #set device
    device = 'cpu'

    # load the data converter which will convert the data to machine-readable format for the algorithm
    converter = spk.data.AtomsConverter(device=device)


    data = open(save_filepath,mode='a',encoding='utf-8-sig')
    save_filepathae= save_filepath.replace('.csv','ae.csv')
    dataae = open(save_filepathae,mode='a',encoding='utf-8-sig')

#    label_n = genfromtxt(labelnfilepath,delimiter=',')
    label_o = genfromtxt(labelofilepath,delimiter=',')
    pk1_data = genfromtxt(pkafilepath,delimiter=',')

    countlabelo = 0
    countlabeln = 0
    for idx in range(start,end):
        if idx%500 == 0:
            print(idx)

        #get molecule properties
        at, props = dataset.get_properties(idx)

        #convert to machine readable format
        inputs = converter(at)

        #use model on inputs  
        pred = model(inputs)

        #get atomic goordinates
        positions = props['_positions']

#        print('idx',idx)
#        utils_xyzdb.usedb(data_filepath,idx,qm9=False,available_properties=available_properties)
        

        from schnetpack.representation.schnet import xs
        from schnetpack.atomistic.output_modules import yi


        xemb = xs[layer].cpu().detach().numpy()
        yi  = yi.detach().numpy()

        number_atoms = len(positions)
        atomic_numbers = props['_atomic_numbers']

        allowed_o_labels = {0,3,4,5,7,8,9,10,11,12,14,16,17,21,26,27,28,29}

        for i in range(number_atoms):
#            if props['_atomic_numbers'][i] == 7:
                
#                if label_n[countlabeln][0] == 28:
#                    print(idx)

                #if label_n is NO2 skip
#                if label_n[countlabeln][0] in allowed_n_labels:

#                    row = ''
#                    for k in range(n_features):
#                        row = row + str(xemb[0][i][k]) + ', '
#                    row = row + '\n'

#                    rowae = '' + str(pk1_data[idx][1]) +', ' + str(label_n[countlabeln][0]) + '\n'

#                    print('pkan',rowae)

#                    data.write(row)
#                    dataae.write(rowae)
#                countlabeln = countlabeln + 1

            if props['_atomic_numbers'][i] == 8:

 #               print('labelo',label_o[countlabelo][0])
                #if label_o is any aldehyde, skip
                if label_o[countlabelo][0] in allowed_o_labels:

                    row = ''
                    for k in range(n_features):
                        row = row + str(xemb[0][i][k]) + ', '
                    row = row + '\n'

                    rowae = '' + str(pk1_data[idx][1]) +', ' + str(label_o[countlabelo][0]) + '\n'

#                    print('pkao',rowae)

                    data.write(row)
                    dataae.write(rowae)       

                countlabelo = countlabelo + 1


    data.close()