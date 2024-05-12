import numpy as np
from numpy import genfromtxt, savetxt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import  PCA
from sklearn.manifold import TSNE

import scipy.linalg as la
import pandas as pd
'''
Last Updated: 2024-03-09

Dimension reduction tools
    1) pca      - an sklearn pca wrapper function to handle pca fitting and testing 
                  straight from datafiles with control options prespecified
    2) tsne     - an sklearn tsne wrapper function to hande tsne fitting straight from datafiles
'''

def pca(FITDATA_FILEPATH,N_COMPONENTS,N_FEATURES,SCALE_DATA,SAVE_FILEPATH,HEADER_EXISTS_FIT,HEADER_EXISTS_APPLY,APPLYDATA_FILEPATH):
    '''
    pca fitting and testing, can apply pca fit on a different set than the one used for fitting

        Args:
            FITDATA_FILEPATH       - the csv dataset that will be used for making the pca fit
            APPLYDATA_FILEPATH     - the csv dataset that will be used for applying the pca fit 
                                     (can be the same as the one used for the fit)
            N_COMPONENTS           - number of components for the pca fit
            N_FEATURES             - number of features in the data
            SCALE_DATA             - boolean to decide to scale the data (StandardScalar) 
                                     before fitting/applying
            SAVE_FILEPATH          - where to save the results of the pca on the APPLY data
            HEADER_EXISTS_FIT      - 1/0 for np.genfromtxt to skip header on the fit
                                     if there is a header then 1
            HEADER_EXISTS_APPLY    - 1/0 for np.genfromtxt to skip header on the apply


        Process:
            fit_data                - loading the data as np.array that will be used for making the pca fit
            apply_data              - loading the data as np.array that will be used for applying the fit
            n_data                  - number of data points for the fitting
            x_data                  - narrows data down to the features of fit data because it also contains labels,
                                      so only N_FEATURES should be permitted
            scaler                  - initiating the StandardScalar, if SCALE_DATA = True, will be used to scale data
                                      based on its meand and variance
            pca                     - initiating the pca fitting class, a constant random seed, and specifying the components 
                                      of the fit
            pca.fit                 - fitting the pca on the data (pca.fit)
            x_pca                   - the applied pca results on the apply data (pca.transform)
            cov                     - the covariance matrix used in the pca fit (measures the covariance between features in the data, 128x128 matrix)
            eig                     - the eigenvalue spectrum of the PCA fit, we will take the log of this for visualization on a plot
            ev                      - the eigenvectors of the PCA fit                
            save_filepathpca        - the filepath of the PCA will be the same as the data but replaced with 'pca.csv'
                                      it will be labelled according to data's label (hstacked right at the end of the PCA projection) 
            save_filepatheig        - filepath of the log(eig) will be the same as the data but replaced with 'eig.csv'
            save_filepathcov        - filepath of the cov will be the same as the data but replaced with 'cov.csv'
            save_filepathev         - filepath of the ev will be the same as the data but replaced with 'ev.csv'

        Returns:
            the pca projection on the apply data, saved to a file, the log eigenvalues of the fit,
            the eigenvectors, and the covariance matrix all saved in files labelled as such
    '''
    
    if HEADER_EXISTS_FIT == False:
        fit_data = pd.read_csv(FITDATA_FILEPATH,delimiter=',',index_col=False,header=None)
    else:
        fit_data = pd.read_csv(FITDATA_FILEPATH,delimiter=',',index_col=False)

    if HEADER_EXISTS_APPLY == False:
        apply_data = pd.read_csv(APPLYDATA_FILEPATH,delimiter=',',index_col=False,header=None)
    else:
        apply_data = pd.read_csv(APPLYDATA_FILEPATH,delimiter=',',index_col=False)

    n_data = len(fit_data)
    x_data = fit_data.iloc[0:n_data,0:N_FEATURES].values

    print(x_data.shape)
    if SCALE_DATA == True:
        scaler = StandardScaler()
        scaler.fit(x_data)
        x_data = scaler.transform(x_data,random_state=100)

    #perform PCA decomposition of the data
    pca = PCA(random_state=100,n_components=N_COMPONENTS)
    pca.fit(x_data)


    x_pca = pca.transform(apply_data.iloc[:,0:N_FEATURES].values)

    cov = pca.get_covariance()
    eig, ev = la.eig(cov)
#    total = sum(eig.real)
#    evt=np.transpose(ev)
#    unit=np.matmul(ev,evt)

    x_pca = pd.DataFrame(x_pca)

    #add columns after N_FEATURES 
    x_pca = pd.concat([x_pca,apply_data.iloc[:n_data,N_FEATURES:]],axis=1)
    eig = np.column_stack((np.arange(0,N_FEATURES,1),np.log(eig),eig))
    
    
    save_filepathpca = SAVE_FILEPATH.replace('.csv','pca.csv')
    x_pca.to_csv(save_filepathpca,index=False)

    save_filepatheig = SAVE_FILEPATH.replace('.csv','eig.csv')
    eig=pd.DataFrame(eig.real)
    eig.to_csv(save_filepatheig,index=False)


    save_filepathev = SAVE_FILEPATH.replace('.csv','ev.csv')
    ev=pd.DataFrame(ev.real)
    ev.to_csv(save_filepathev,index=False)

    save_filepathcov = SAVE_FILEPATH.replace('.csv','cov.csv')
    cov=pd.DataFrame(cov.real)
    cov.to_csv(save_filepathcov,index=False)




def tsne(DATA_FILEPATH, DIMENSION, PERP, N_FEATURES, N_DATA):
    '''
    tsne fitting for visualization of high dimensional data in low dimensions

        Args:
            DATA_FILEPATH       - the filepath of the high dimensional data
            DIMENSION           - the dimension of the projection, 2D or 3D
            PERP                - the perplexity of the projection 
                                  (min number of effective neighbors per nbrhood)
            N_FEATURES          - number of features in the data 
            N_DATA              - number of data in the filepath
        Process:
            full_data           - the full data set with its labels
            x_data              - only the features part of the data (not the label)
            X                   - tsne projection using dimension and perplexity on x_data
            save_name           - tsne will be saved with the perplexity number as 'tsne#perp.csv'
                                  to keep track of perplexity 
            save_filepath       - save filepath will be the same as the data
                                  but replaced with 'tsne#perp.csv'

        Returns
            saved csv of the tsne projection, 
            labelled according to data's label
    '''

    full_data = genfromtxt(DATA_FILEPATH,delimiter=',',encoding='utf-8-sig',skip_header=1)
    x_data = full_data[0:N_DATA,0:N_FEATURES]

    

    X = TSNE(n_components=DIMENSION,perplexity=PERP).fit_transform(x_data)

    X = np.hstack((X,full_data[:,N_FEATURES:]))


    save_name = 'tsne%s.csv' %(PERP)
    save_filepath = DATA_FILEPATH.replace('.csv', save_name)
    print(save_filepath)
    savetxt(save_filepath,X,delimiter=',')

