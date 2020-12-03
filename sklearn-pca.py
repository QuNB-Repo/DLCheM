# The input of this program is a .csv data sheet
# set up in a special format (see given example). 
# The program standardizes the input data then 
# performs PCA analysis on it. Outputs PCA graph 
# and covariance matrix. 

#import libraries

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# read excel (examples of how to set up your data given in folder)
data = pd.read_csv("./PCACME-9.csv") 

# scale the data 
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
scaler.fit(data)
scaled_data = scaler.transform(data)

#print(scaled_data)

#perform PCA fit and decompose the data
from sklearn.decomposition import  PCA
pca = PCA(n_components = 10)
pca.fit(scaled_data)
x_pca = pca.transform(scaled_data)      

#print(x_pca.shape)

#label= pd.read_csv("./label.csv")

#plot figure
plt.figure(figsize=(8,6))
plt.scatter(x_pca[:,3],x_pca[:,5])
plt.xlabel('First principal component')
plt.ylabel('Second principal component')


cov = pca.get_covariance()
#print(cov)
# optional print covariance and associated eigenvectors
#print('covariance')
#print(cov)

#print(pca.components_) 

import scipy.linalg as la
eig, ev = la.eig(cov)
print(eig)
#print('all eigenvalues')
#print(eig) 
#print(ev[:,0])
