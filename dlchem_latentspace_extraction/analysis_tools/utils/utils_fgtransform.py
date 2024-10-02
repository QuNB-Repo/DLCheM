import numpy as np
from numpy import savetxt, genfromtxt
import fileinput

def labelfile(n_data,label_id,save_filepath):

    label = np.zeros((n_data,2))
    idx = 0
    for i in range(n_data):
        label[i][0] = label_id
        label[i][1] = idx  
        idx = idx + 1

    savetxt(save_filepath,label,delimiter=',')


def vecdiffmean(data_filepath1,data_filepath2,num_features,n_molecules,fg_trans_fildir):

    data1 = genfromtxt(data_filepath1,delimiter=',')
    data2 = genfromtxt(data_filepath2,delimiter=',')
 
    vec_diff_matrix = np.zeros((n_molecules,num_features))
    for i in range(n_molecules):
        for j in range(num_features):
            vec_diff_matrix[i][j] = data1[i][j] - data2[i][j]
    
    vec_diff_mean = np.mean(vec_diff_matrix,axis=0)

    save_path = fg_trans_fildir + 'diff.csv'
    savetxt(save_path,vec_diff_mean)
    return vec_diff_mean

def vecdiffmean_nonlinear(data_filepath1,data_filepath2,num_features,n_molecules,fg_trans_fildir):

    data1 = genfromtxt(data_filepath1,delimiter=',')
    data2 = genfromtxt(data_filepath2,delimiter=',')
 
    vec_diff_matrix = np.zeros((n_molecules,num_features))
    for i in range(n_molecules):
        for j in range(num_features):
            vec_diff_matrix[i][j] = data1[i][j] - np.tanh(data2[i][j])

    vec_diff_mean = np.mean(vec_diff_matrix,axis=0)

    save_path = fg_trans_fildir + 'Odiff.csv'
    savetxt(save_path,vec_diff_mean)
    return vec_diff_mean



def add_vectomat(vec_filepath,data_filepath,n_features,save_filepath):

    vec = genfromtxt(vec_filepath,delimiter=',',encoding='utf-8-sig')
    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')
    print('data',data)
    print('vec',vec)

    new_data = np.zeros((len(data),len(vec)))
    for i in range(len(data)):
        for j in range(n_features):
            new_data[i][j] = data[i][j] - vec[j]

    new_data = np.hstack((new_data,data[:,n_features:]))

    savetxt(save_filepath,new_data,delimiter=',')

def add_vectomat_tanh(vec_filepath,data_filepath,n_features,save_filepath):

    vec = genfromtxt(vec_filepath,delimiter=',',encoding='utf-8-sig')
    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')
    print('data',data)
    print('vec',vec)

    new_data = np.zeros((len(data),len(vec)))
    for i in range(len(data)):
        for j in range(n_features):
            new_data[i][j] = data[i][j] - vec[j]

    new_data = np.hstack((new_data,data[:,n_features:]))

    savetxt(save_filepath,new_data,delimiter=',')



def nrst_vec(vec_filepath,data_filepath):
    vec = genfromtxt(vec_filepath,delimiter=',',encoding='utf-8-sig')
    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')

    pred_ind = np.zeros(len(vec))
    for i in range(len(vec)):
        nrst_dist = 100
        dist = 0
        for j in range(len(data)):
            dist = np.linalg.norm(data[j,:]-vec[i,:])
            if dist < nrst_dist:
                nrst_dist = dist
                pred_ind[i] = j
    
    print('pred_ind \n',pred_ind)
    error = 0
    for i in range(len(pred_ind)):
        integer = int(pred_ind[i])
        if integer != i:
            error = error + 1
    
    error = (error/len(pred_ind))*100
    print('error \n', error)

def mean_vec(data_filepath):
    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')
    mean_vec = np.mean(data,axis=0)
    mean_size = np.linalg.norm(mean_vec)

    print(mean_vec)
    print(mean_size)
    return mean_vec, mean_size

def nearest_trans(true_trans_filepath,art_trans_filepath,n_molecules,n_features,start_trans_idx=0,stack_qm9=False,qm9embs_filepath=''):
    true_trans = genfromtxt(true_trans_filepath,delimiter=',')
    art_trans = genfromtxt(art_trans_filepath,delimiter=',')
    
    if stack_qm9 == True:
        qm9_embs = genfromtxt(qm9embs_filepath,delimiter=',')
        true_trans = np.vstack((true_trans,qm9_embs))
        n_molecules_full = len(true_trans)
    else:
        n_molecules_full = len(true_trans)


    neighbor = np.zeros(n_molecules)
    match = 0
    for i in range(n_molecules):
        if i % 1000 == 0:
            print(i)

        nearest = 1000

        for j in range(0,n_molecules):
            distance = np.linalg.norm(art_trans[i,0:n_features] - true_trans[start_trans_idx+j,0:n_features])
            if distance < nearest:
                nearest = distance 
                neighbor[i] = j
        if neighbor[i] == i:
            match = match + 1
        elif neighbor[i] != i:
            print(i)
    
    savetxt(art_trans_filepath.replace('.csv','nbrs.csv'),neighbor,delimiter=',')
    print('matched',(match/n_molecules)*100,'%')
    


def scatters_to_vectors(data_filepath,save_filepath):

    #generate data
    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')

    vector_data = np.zeros((int(len(data)/2),4))
    #Assuming hald the data is init, half is final
    for each_data in range(0,int(len(data)/2)):
        vector_data[each_data,0:2] = data[each_data,0:2]
        vector_data[each_data,2:4] = data[each_data+int(len(data)/2),0:2]

    print(vector_data)
    savetxt(save_filepath,vector_data,delimiter=',')


