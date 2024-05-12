import numpy as np
from numpy import genfromtxt, savetxt

import matplotlib.pyplot as plt

from schnetpack.datasets import QM9

import pandas as pd

def euc_dist_bw_layers(layer1_filepath,layer2_filepath):

    #generate data from files
    layer1 = pd.read_csv(layer1_filepath,delimiter=',')
    layer2 = pd.read_csv(layer2_filepath,delimiter=',')

    #generate data from files
#    layer1 = pd.read_csv(layer1_filepath,delimiter=',')
#    layer2 = pd.read_csv(layer2_filepath,delimiter=',')

    fg_mask = layer1['ldalabel'] == 4
    fg_mask2 = layer1['ldalabel'] == 5

    embs1 = layer1.loc[fg_mask, 'embs0':'embs127']
    embs2 = layer1.loc[fg_mask2, 'embs0':'embs127']

    print('HERE',embs1.iloc[0])
#    embs1 = layer1.iloc[0:,0:128]
#    embs2 = layer2.iloc[0:,0:128]

    distance = [np.linalg.norm(embs1.iloc[each_datapoint]-embs2.iloc[each_datapoint]) for each_datapoint in range(len(embs2))]

    average_distance = np.mean(distance)

    return distance, average_distance




def distmatXvsYscatter(data_filepath,save_filepath,x_features,y_features):

    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')

    out_data = []
    for datapoint1 in range(len(data)):
        target = data[datapoint1,x_features[0]:x_features[1]]

        for datapoint2 in range(len(data)):
        
            compare_with = data[datapoint2,x_features[0]:x_features[1]]

            distance = np.linalg.norm(target-compare_with)

#        labelpoint = label[datapoint][0]

            out_data.append([distance])
    
    y_data = []
    for datapoint1 in range(len(data)):
        target = data[datapoint1,y_features[0]]

        y_of_this = data[datapoint1][y_features[0]]
        print('here',y_of_this)

        for datapoint2 in range(len(data)):
        
            compare_with = data[datapoint2,y_features[0]]

            compare_y = data[datapoint2][129]
            print('compare',compare_y)
            if y_of_this == compare_y:
                y_label = '0'
            else:
                y_label = '12632256'




            distance = np.linalg.norm(target-compare_with)

#        labelpoint = label[datapoint][0]
            y_data.append([distance,y_label])
    
    out_data = np.asarray(out_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    out_data = np.column_stack((out_data,y_data))
    
    print(out_data)
    savetxt(save_filepath,out_data,delimiter=',',encoding='utf-8-sig')


def euc_dist_vs_atomy(element_idx,data_filepath,save_filepath):

    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')

    target = data[element_idx,0:128]

    out_data = []
    for datapoint in range(len(data)):
        
        compare_with = data[datapoint,0:128]

        distance = np.linalg.norm(target-compare_with)

#        labelpoint = label[datapoint][0]

        out_data.append([distance])
    
    savetxt(save_filepath,out_data,delimiter=',')




def EucDistDistribution(element_idx,label_filepath,data_filepath,bin_size,numb_bins,init_bin,numb_categories,save_filepath):

    label = genfromtxt(label_filepath,delimiter=',')
    data = genfromtxt(data_filepath,delimiter=',')

    #for each datapoint, for each category, calculate distance to element_idx in that category


    colors=['red','orange','yellow','mediumaquamarine','mediumseagreen',
            'palegreen','crimson','lightblue','cadetblue','lightskyblue',
            'royalblue','cornflowerblue','khaki','gold','firebrick',
            'lightpink','mistyrose','deeppink','hotpink','mediumvioletred',
            'plum','purple','palevioletred','mediumorchid','brown',]
    colors=['red','orange','yellow','mediumaquamarine','mediumseagreen',
            'palegreen','yellowgreen','lightblue','cadetblue','lightskyblue',
            'royalblue','cornflowerblue','khaki','gold','peru',
            'lightpink','mistyrose','deeppink','hotpink','mediumvioletred',
            'plum','purple','palevioletred','mediumorchid','tan',
            'indianred','lightsalmon','lightcoral']
    all_dists = []

    
    for category in range(numb_categories[0],numb_categories[1]):
        euc_dist = []
        for datapoint in range(len(data)):
            if label[datapoint][0] == category: 
                euc_dist.append([np.linalg.norm(data[datapoint,:]-data[element_idx,:]),datapoint,label[datapoint][1]])
    
    
        bin = init_bin
        dist = np.zeros((numb_bins,2))
        for i in range(numb_bins):
            dist[i][0] = bin 
            for j in range(len(euc_dist)):
                if bin <= euc_dist[j][0] < bin + bin_size:
                    dist[i][1] = dist[i][1] + 1
            bin = bin + bin_size

        all_dists.append(dist[:,1])
        plt.plot(dist[:,0],dist[:,1],color=colors[category])
    plt.show()

    all_dists.append(dist[:,0])
    all_dists = np.transpose(all_dists)
    np.savetxt(save_filepath,all_dists,delimiter=',')

    return euc_dist


def extractrange(euc_dist,label,min,max):
    
    element_idx = []
    for i in range(len(euc_dist)):
        if min <= euc_dist[i][0] < max:
            print('yessss')
            element_idx.append(euc_dist[i][1])
    
    return element_idx


        

def distvsE(data_filepath,label_filepath,qm9_filepath,element_idx):

    data = genfromtxt(data_filepath,delimiter=',')
    label = genfromtxt(label_filepath,delimiter=',')

    qm9data = QM9(qm9_filepath,download=False,remove_uncharacterized=True)
    ref_idx = int(label[element_idx][1])
    at, props = qm9data.get_properties(ref_idx)

    ref_datapoint = data[element_idx,:]
    ref_energy =  props[QM9.U0][0]

    dists = []
    for index, datapoint in enumerate(data):

        dist = np.linalg.norm(datapoint-ref_datapoint)

        molecule_idx = int(label[index][1])

        at, props = qm9data.get_properties(molecule_idx)
        energy =  props[QM9.U0][0]

        diff_energy = abs(ref_energy - energy)

        dists.append([dist,diff_energy])
    
    return dists


def distvsAE(data_filepath,ae_filepath,element_idx):

    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')
    ae_data = genfromtxt(ae_filepath,delimiter=',',encoding='utf-8-sig')

    ref_datapoint = data[element_idx,:]
    ref_ae =  ae_data[element_idx]

    dists = []
    for index, datapoint in enumerate(data):

        dist = np.linalg.norm(datapoint-ref_datapoint)

        ae =  ae_data[index]

        diff_ae = abs(ref_ae - ae)

        dists.append([dist,diff_ae])
    
    return dists



