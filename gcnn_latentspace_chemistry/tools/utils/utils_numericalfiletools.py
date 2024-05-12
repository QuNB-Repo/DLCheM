import numpy as np
from numpy import savetxt, genfromtxt

def vstacktwofiles(filepath1,filepath2,save_filepath,n_data1,n_data2,single_dim,skip_header1,skip_header2):
    data1 = np.genfromtxt(filepath1,delimiter=',',encoding='utf-8-sig',skip_header=skip_header1)
    data2 = np.genfromtxt(filepath2,delimiter=',',encoding='utf-8-sig',skip_header=skip_header2)
    if single_dim == True:
        data = np.concatenate([data1[:n_data1,:],data2[:n_data2,:]])
    else: 
        data = np.vstack((data1,data2))

    savetxt(save_filepath,data,delimiter=',')

def hstacktwofiles(filepath1,filepath2,single_dim,save_filepath):
    data1 = np.genfromtxt(filepath1,delimiter=',',encoding='utf-8-sig')
    data2 = np.genfromtxt(filepath2,delimiter=',',encoding='utf-8-sig')
    if single_dim == True:
        data = np.concatenate([data1,data2])
    else: 
        data = np.hstack((data1,data2))

    savetxt(save_filepath,data,delimiter=',')

def removelastcol(filepath):

    data = np.genfromtxt(filepath,delimiter=',',encoding='utf-8-sig')

    newdata = data[:,0:len(data)-1]

    savetxt(filepath,newdata,delimiter=',')



def split_2(dataset_filepath,save_filepath):

    data = np.genfromtxt(dataset_filepath,delimiter=',',encoding='utf-8-sig')

    data_first = data[0:int(len(data)/2)]

    data_second = data[int(len(data)/2):int(len(data))]
    
    save_filepath_first = save_filepath + '1.csv'
    save_filepath_second = save_filepath + '2.csv'

    savetxt(save_filepath_first,data_first,delimiter=',')
    savetxt(save_filepath_second,data_second,delimiter=',')


def addtwofiles(dataset1_filepath,datset2_filepath,features,save_filepath):
    # Load the first data file
    data1 = np.genfromtxt(dataset1_filepath, skip_header=1, delimiter=',', usecols=range(features[0], features[1]))

    # Load the second data file
    data2 = np.genfromtxt(datset2_filepath, skip_header=1, delimiter=',', usecols=range(features[0], features[1]))

    # Add the two arrays element-wise
    sum_data = data1 + data2

    # Print the resulting array
    savetxt(save_filepath,sum_data,delimiter=',')