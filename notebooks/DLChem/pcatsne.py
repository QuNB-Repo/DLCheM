# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 18:02:03 2021

@author: aelsamma
"""

import numpy as np
from numpy import genfromtxt, savetxt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import  PCA
from sklearn.manifold import TSNE

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d

import pandas as pd

import scipy.linalg as la


from pylab import *
import matplotlib.pyplot as plt
import matplotlib.lines
from matplotlib.transforms import Bbox, TransformedBbox
from matplotlib.legend_handler import HandlerBase
from matplotlib.image import BboxImage


class HandlerLineImage(HandlerBase):

    def __init__(self, path, space=0, offset = 0 ):
        self.space=space
        self.offset=offset
        self.image_data = plt.imread(path)        
        super(HandlerLineImage, self).__init__()

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):

        l = matplotlib.lines.Line2D([xdescent+self.offset,xdescent+(width-self.space)/3.+self.offset],
                                     [ydescent+height/2., ydescent+height/2.])
        l.update_from(orig_handle)
        l.set_clip_on(False)
        l.set_transform(trans)

        bb = Bbox.from_bounds(xdescent +(width+self.space)/3.+self.offset,
                              ydescent,
                              height*self.image_data.shape[1]/self.image_data.shape[0],
                              height)

        tbb = TransformedBbox(bb, trans)
        image = BboxImage(tbb)
        image.set_data(self.image_data)

        self.update_prop(image, orig_handle, legend)
        return [l,image]


def pca(element,file_path,label_path,scale_data,n_components,dimension):
    if element == 'H':
         colors=['purple','none','palegreen','pink','lightblue',
             'mediumseagreen','lightskyblue','lightcoral','none','darkseagreen',
             'cornflowerblue','indianred','none','yellowgreen','palevioletred',
             'aquamarine','none','steelblue','mediumvioletred','none',
             'darkolivegreen','none','none','khaki','hotpink',
             'navajowhite','none','none','thistle','darksalmon',
             'tomato','none','burlywood']
    if element == 'C':
        colors=['black','none','#ccffcc','#ccffff','#ffcccc',
             '#99ff66','none','none','#66ffff','#ff9999',
             '#77b300','#33cccc','#ff6666','none','#666633',
             '#cc0000','none','none','none','none',
             'none','#ffffcc','#ffcc99','#66ccff','#ccccff',
             'none','#00cc66','#009999','#ff9933','#0099ff',
             '#9999ff','black','#00ff99','#006666','black',
             'none','#99ccff', '#3333cc', '#ff6600','none',
             'blue','none','none','black','black',
             'black','none']
    if element == 'N':
        colors=['none','none','none','palegreen','khaki',
             'none','mediumaquamarine','mediumseagreen','none','orange',
             'darkseagreen','none','none','none','none',
             'lightskyblue', 'steelblue','lightcoral','none','none',
             'tomato','plum','purple','hotpink']
    
    if element == 'O':
        colors=['none','none','tomato','none','mediumaquamarine',
             'mediumseagreen','palegreen','none','lightblue','none',
             'lightskyblue','none','cornflowerblue','khaki','NONE',
             'none','lightpink','mistyrose','deeppink','hotpink',
             'mediumvioletred','plum','purple','NONE','mediumorchid']

    data = genfromtxt(file_path,delimiter=',')
    print(data)
    
    if scale_data == True:
        scaler = StandardScaler()
        scaler.fit(data)
        data = scaler.transform(data)
    #perform PCA decomposition of the data

    
    pca = PCA(n_components)
    pca.fit(data)
    x_pca = pca.transform(data)
    
    cov = pca.get_covariance()
    eig, ev = la.eig(cov)
    total = sum(eig.real)
    #plot the PCA if dimension 2 or 3, if element H or O for colors/legend
    x=x_pca
    pca_file = '../../../data/pca/%s.csv' %(str(number_data))
    savetxt(pca_file, x, delimiter=',')
    evt=np.transpose(ev)
    unit=np.matmul(ev,evt)

    
    label = pd.read_csv(label_path,delimiter=',')

    fig4,ax4 = plt.subplots(figsize=(6,10),dpi=1000,nrows=2,ncols=1)
    
    ax4[1].set_xlabel('t-SNE1')
    ax4[1].set_ylabel('t-SNE2')
    

    perp = 50
    dimension=2
    
    print('here')
    X = TSNE(n_components=dimension,perplexity=perp).fit_transform(data)
                 
    
    if dimension == 2:
        ax4[1].scatter(X[:,0],X[:,1],c=label['Target'],cmap=matplotlib.colors.ListedColormap(colors))
    if dimension == 3:
        ax4[1].scatter(X[:,0],X[:,1],X[:,2],c=label['Target'],cmap=matplotlib.colors.ListedColormap(colors))
    
    ax4[0].set_xlabel('First PC')
    ax4[0].set_ylabel('Second PC')
    ax4[0].scatter(x[:,0],x[:,1],c=label['Target'],cmap=matplotlib.colors.ListedColormap(colors))
        
        
    line1,  = plt.plot(0,0, color='tomato', lw=10)
    line2,  = plt.plot(0,0, color='mediumaquamarine', lw=10)  
    line3, = plt.plot(0,0, color='mediumseagreen', lw=10)  
    line4, = plt.plot(0,0, color='palegreen', lw=10)  
    line5, = plt.plot(0,0, color='lightblue', lw=10)  
    line6, = plt.plot(0,0, color='lightskyblue', lw=10)  
    line7, = plt.plot(0,0, color='cornflowerblue', lw=10)  
    line8, = plt.plot(0,0, color='khaki', lw=10)  
    line9, = plt.plot(0,0, color='lightpink', lw=10)  
    line10, = plt.plot(0,0, color='lightcoral', lw=10)  
    line11, = plt.plot(0,0, color='deeppink', lw=10)  
    line12, = plt.plot(0,0, color='hotpink', lw=10)  
    line13, = plt.plot(0,0, color='mediumvioletred', lw=10)  
    line14, = plt.plot(0,0, color='plum', lw=10)  
    line15, = plt.plot(0,0, color='purple', lw=10)  
    line16, = plt.plot(0,0, color='mediumorchid', lw=10)  
    
    
    
    if element == 'O':
        l1 = ax4[0].legend([line1,line2,line3,line4,line5,line6,line7,line8], ["","","","","","","",""],
           handler_map={ line1: HandlerLineImage(image_dir + "./oxy_COC.png"), line2: HandlerLineImage(image_dir +"./oxy_HOC-CCC.png"),
                        line3: HandlerLineImage(image_dir + "./oxy_HOC-CCH.png"), line4: HandlerLineImage(image_dir +"./oxy_HOC-CHH.png"),
                        line5: HandlerLineImage(image_dir + "./oxy_HOC-CC.png"),   line6: HandlerLineImage(image_dir +"./oxy_HOC-CN.png"),
                        line7: HandlerLineImage(image_dir + "./oxy_HOC-NN.png"), line8: HandlerLineImage(image_dir +"./oxy_HON.png")},
           handlelength=2, labelspacing=0.0, fontsize=36, borderpad=0.15, loc=3, 
            handletextpad=0.2, borderaxespad=0.15,bbox_to_anchor=(1.0, -1.5),prop={'size': 80})
        
        l2 = ax4[1].legend([line9,line10,line11,line12,line13,line14,line15,line16], ["","","","","","","",""],
           handler_map={line9: HandlerLineImage(image_dir + "./oxy_OC-CC.png"), line10: HandlerLineImage(image_dir + "./oxy_OC-CH.png"), 
                        line11: HandlerLineImage(image_dir + './oxy_OC-CN.png'), line12: HandlerLineImage(image_dir + "./oxy_OC-CO.png"), 
                        line13: HandlerLineImage(image_dir + "./oxy_OC-HN.png"), line14: HandlerLineImage(image_dir + "./oxy_OC-NO.png"), 
                        line15: HandlerLineImage(image_dir + "./oxy_OC-NN.png"), line16: HandlerLineImage(image_dir + "./oxy_OC-HO.png")},
           handlelength=2, labelspacing=0.0, fontsize=36, borderpad=0.15, loc=3, 
            handletextpad=0.2, borderaxespad=0.15,bbox_to_anchor=(1.6, -0.3),prop={'size': 80})
    if element == 'N':
        l1 = ax4[0].legend([line1,line2,line3,line4,line5,line6], ["","","","","",""],
           handler_map={ line1: HandlerLineImage(image_dir + "./nitro_NH2-C4.png"), line2: HandlerLineImage(image_dir + "./nitro_NH2-C3.png"),
                        line3: HandlerLineImage(image_dir + "./nitro_NH-C4.png"), line4: HandlerLineImage(image_dir + "./nitro_NH-C3.png"),
                        line5: HandlerLineImage(image_dir + "./nitro_NH-CN.png"),   line6: HandlerLineImage(image_dir + "./nitro_NH-CCC.png")},
           handlelength=2, labelspacing=0.0, fontsize=36, borderpad=0.15, loc=3, 
            handletextpad=0.2, borderaxespad=0.15,bbox_to_anchor=(1.0, -1.2),prop={'size': 80})
        
        
        l1 = ax4[1].legend([line7,line8,line9,line10,line11,line12], ["","","","","",""],
           handler_map={line7: HandlerLineImage(image_dir + "./nitro_N-CC.png"), line8: HandlerLineImage(image_dir + "./nitro_N-HC.png"),
                        line9: HandlerLineImage(image_dir + "./nitro_N-CO.png"), line10: HandlerLineImage(image_dir + "./nitro_N-NN.png"),
                        line11: HandlerLineImage(image_dir + "./nitro_N-C.png"),line12: HandlerLineImage(image_dir+"./nitro_NH2+-C2.png")},
            handlelength=2, labelspacing=0.0, fontsize=36, borderpad=0.15, loc=3, 
            handletextpad=0.2, borderaxespad=0.15,bbox_to_anchor=(1.6, 0.0),prop={'size': 80})

        
    plt.show()
     

element = 'N'
name_data='rep-0'
number_data=5000
file_path = '../../../data/schnet/5000-mod/' + name_data + '%s.csv' %(element)
label_path = '../../../data/label_qm9/%s/5000/label%s.csv' %(element,element) 
scale_data = True
n_components=30
dimension=2
image_dir = '../../../data/label_qm9/%s/' %(element)
pca(element,file_path,label_path,scale_data,n_components,dimension)

