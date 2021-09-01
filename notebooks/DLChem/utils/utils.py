# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 19:45:49 2021

@author: aelsamma
"""

import os

import sys

import fileinput


from pylab import *
import matplotlib.pyplot as plt
import matplotlib.lines
from matplotlib.transforms import Bbox, TransformedBbox
from matplotlib.legend_handler import HandlerBase
from matplotlib.image import BboxImage

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d

import pandas as pd

from numpy import genfromtxt, savetxt
import numpy as np

#Writes xyz file from a loaded database (db), needs props of db as input, and idx
#required for labelling code
def write_xyz_from_db(props,idx,file_path):
    
    name_xyz = file_path+str(idx)+'.xyz'
    # open an empty file
    xyz_file = open(name_xyz,mode='w',encoding='utf-8')

    # copy props['_atomic_numbers'], props['_positions] tensor and change tensor to numpy array 
    atomic_numbers = props['_atomic_numbers']
    atomic_numbers = atomic_numbers.detach().numpy()
    
    number_atoms = len(atomic_numbers)

    positions = props['_positions']
    positions = positions.detach().numpy()
    # write xyz file in xyz file format
    xyz_file.write(str(number_atoms)+'\n')
    xyz_file.write('Title'+'\n')
    for i in range(number_atoms):
        if atomic_numbers[i] == 1:
            xyz_file.write('H ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 6:
            xyz_file.write('C ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 7:
            xyz_file.write('N ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 8:
            xyz_file.write('O ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
        if atomic_numbers[i] == 9:
            xyz_file.write('F ' + str(positions[i][0]) + ' ' + str(positions[i][1]) + ' ' + str(positions[i][2]) + '\n')
    xyz_file.close()
    return name_xyz

#runs obabel conversion tool on the created xyz file to convert to mol file
#required for labelling code
def xyz_to_mol(idx,file_path):
    # define name of temporary xyz file according to idx, define the name of temporary mol file 
    name_xyz = file_path + str(idx) + '.xyz'
    name_mol = file_path + str(idx) + '.mol'
    
    #use obabel to convert xyz to mol
    os.system('obabel ' + name_xyz + ' -O ' + name_mol)
    
    return name_mol

def get_platform():
    platforms = {
        'linux1' : 'Linux',
        'linux2' : 'Linux',
        'darwin' : 'OS X',
        'win32' : 'Windows'
    }
    if sys.platform not in platforms:
        return sys.platform
    
    return platforms[sys.platform]



def store_positions(idx,name_mol,props,element):
    countmollines = 0
    array_of_element = []
    for line in fileinput.FileInput(name_mol,inplace=0):
        countmollines = countmollines + 1
        if element in line:
            array_of_element.append([])
            array_of_element[len(array_of_element)-1] = int(countmollines-4)
    return array_of_element


def connection_matrix(array_of_element,name_mol,number_atoms):
    count=0
    connected_positions = []
    array = []
    for j in range(len(array_of_element)):
        connected_positions.append([])
        array.append([])
        array[j] = str(array_of_element[j])
    for line in fileinput.FileInput(name_mol,inplace=0):
        count = count + 1
        if count >= number_atoms + 5:
            line = line.replace(line,line[:7])
            for j in range(len(array_of_element)):
                if ' ' + str(array_of_element[j]) + ' ' in line and 'RAD' not in line:
                    connected_positions[j].append([])
                    if str(line[1]+line[2]+line[3]) == ' '+array[j]+' ' or str(line[0]+line[1]+line[2]+line[3]) == ' '+array[j]+' ':
                        connected = str(line[4:6])
                        connected = int(connected)
                        connected_positions[j][len(connected_positions[j])-1]=connected
                    else:
                        connected = str(line[1]+line[2]+line[3])
                        connected = int(connected)
                        connected_positions[j][len(connected_positions[j])-1]=connected
#    print(connected_positions)
    return connected_positions


def neighboring_connections(name_mol,number_atoms,carbon_position):
    carbon_neighbor_connections = []
    countmollines = 0
    for line in fileinput.FileInput(name_mol,inplace=0):
        countmollines = countmollines + 1
        if countmollines >= number_atoms + 5:
            line = line.replace(line,line[:7])
            if ' ' + str(carbon_position) + ' ' in line and 'RAD' not in line:
                carbon_neighbor_connections.append([])
                if str(line[1]+line[2]+line[3]) == ' '+str(carbon_position)+' ' or str(line[0]+line[1]+line[2]+line[3]) == ' '+str(carbon_position)+' ':
                    connected = str(line[4:6])
                    connected = int(connected)
                    carbon_neighbor_connections[len(carbon_neighbor_connections)-1]=connected
                else:
                    connected = str(line[1]+line[2]+line[3])
                    connected = int(connected)
                    carbon_neighbor_connections[len(carbon_neighbor_connections)-1]=connected
    return carbon_neighbor_connections


def check_O(neighbor,name_xyz):
    countxyzlines = 0
    O_present=False

    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'O ' in line:
                O_present=True
    return O_present


def check_N(neighbor,name_xyz):
    countxyzlines = 0
    N_present=False

    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'N ' in line:
                N_present=True
    return N_present

def check_C(neighbor,name_xyz):
    countxyzlines = 0
    C_present=False

    for line in fileinput.FileInput(name_xyz,inplace=0):
        countxyzlines = countxyzlines + 1
        if countxyzlines == neighbor + 2:
            if 'C ' in line:
                C_present=True
    return C_present

def count_nn(nn,name_xyz):
    countC = 0
    countH = 0
    countO = 0
    countN = 0
    for j in range(len(nn)):
        target = nn[j]
        countxyzlines = 0
        for line2 in fileinput.FileInput(name_xyz,inplace=0):
            countxyzlines = countxyzlines+1
            if countxyzlines == nn[j]+2:
                if 'C ' in line2:
                    countC = countC+1
                if 'H ' in line2:
                    countH = countH+1
                if 'N ' in line2:
                    countN = countN+1
                if 'O ' in line2:
                    countO = countO+1
    total = countC + countH + countO + countN
    
    return total, countC, countH, countO, countN


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
    
    
def plotwlabels(data,label_path,dimension,element,colors,markers,mark,col,image_dir,xlabel,ylabel,legend=False):
    
    data = genfromtxt(data,delimiter=',')
    
    label = pd.read_csv(label_path,delimiter=',')

    fig4,ax4 = plt.subplots(figsize=(10,6),dpi=500,nrows=1,ncols=1)
    
    
    if mark == True:
        for j in range((48)):
            for i in range((len(data))):
    
                if label['Target'][i] == j:
                    if dimension == 2:    
                        ax4.scatter(data[i][0],data[i][1],color=colors[j],marker=markers[j])
   
    elif col == True:
        if dimension == 2:
            ax4.scatter(data[:,0],data[:,1],c=label['Target'],cmap=matplotlib.colors.ListedColormap(colors))
        if dimension == 3:
            ax4 = plt.axes(projection='3d')
            ax4.scatter(data[:,0],data[:,1],data[:,2],c=label['Target'],cmap=matplotlib.colors.ListedColormap(colors))            
#            ax4.view_init(210,320)
    else:
        if dimension == 2:
            ax4.scatter(data[:,0],data[:,1])
            
            
    
    ax4.set_xlabel(xlabel)
    ax4.set_ylabel(ylabel)

#    ax4[0].scatter(data[:,0],data[:,1],c=label['Target'],cmap=matplotlib.colors.ListedColormap(colors))
    
    
    if legend == True:
    
        if element == 'O':
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
            l1 = ax4.legend([line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12,line13,line14,line15,line16], ["","","","","","","","","","","","","","","",""],
               handler_map={ line1: HandlerLineImage(image_dir + "./oxy_COC.png"), line2: HandlerLineImage(image_dir +"./oxy_HOC-CCC.png"),
                            line3: HandlerLineImage(image_dir + "./oxy_HOC-CCH.png"), line4: HandlerLineImage(image_dir +"./oxy_HOC-CHH.png"),
                            line5: HandlerLineImage(image_dir + "./oxy_HOC-CC.png"),   line6: HandlerLineImage(image_dir +"./oxy_HOC-CN.png"),
                            line7: HandlerLineImage(image_dir + "./oxy_HOC-NN.png"), line8: HandlerLineImage(image_dir +"./oxy_HON.png"),
                            line9: HandlerLineImage(image_dir + "./oxy_OC-CC.png"), line10: HandlerLineImage(image_dir + "./oxy_OC-CH.png"), 
                            line11: HandlerLineImage(image_dir + './oxy_OC-CN.png'), line12: HandlerLineImage(image_dir + "./oxy_OC-CO.png"), 
                            line13: HandlerLineImage(image_dir + "./oxy_OC-HN.png"), line14: HandlerLineImage(image_dir + "./oxy_OC-NO.png"), 
                            line15: HandlerLineImage(image_dir + "./oxy_OC-NN.png"), line16: HandlerLineImage(image_dir + "./oxy_OC-HO.png")},
               handlelength=2, labelspacing=0.0, fontsize=36, borderpad=0.15, loc=3, 
                handletextpad=0.2, borderaxespad=0.15,bbox_to_anchor=(1.6, -0.3),prop={'size': 80})
        
        
        if element == 'N':
            line1,  = plt.plot(0,0, color='palegreen', lw=10)
            line2,  = plt.plot(0,0, color='khaki', lw=10)  
            line3, = plt.plot(0,0, color='mediumaquamarine', lw=10)  
            line4, = plt.plot(0,0, color='mediumseagreen', lw=10)  
            line5, = plt.plot(0,0, color='orange', lw=10)  
            line6, = plt.plot(0,0, color='darkseagreen', lw=10)  
            line7, = plt.plot(0,0, color='lightskyblue', lw=10)  
            line8, = plt.plot(0,0, color='steelblue', lw=10)  
            line9, = plt.plot(0,0, color='lightcoral', lw=10)  
            line10, = plt.plot(0,0, color='tomato', lw=10)  
            line11, = plt.plot(0,0, color='plum', lw=10)  
            line12, = plt.plot(0,0, color='hotpink', lw=10)
            l1 = ax4.legend([line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12], ["","","","","","","","","","","",""],
               handler_map={ line1: HandlerLineImage(image_dir + "./nitro_NH2-C4.png"), line2: HandlerLineImage(image_dir + "./nitro_NH2-C3.png"),
                            line3: HandlerLineImage(image_dir + "./nitro_NH-C4.png"), line4: HandlerLineImage(image_dir + "./nitro_NH-C3.png"),
                            line5: HandlerLineImage(image_dir + "./nitro_NH-CN.png"),   line6: HandlerLineImage(image_dir + "./nitro_NH-CCC.png"),
                            line7: HandlerLineImage(image_dir + "./nitro_N-CC.png"), line8: HandlerLineImage(image_dir + "./nitro_N-HC.png"),
                            line9: HandlerLineImage(image_dir + "./nitro_N-CO.png"), line10: HandlerLineImage(image_dir + "./nitro_N-NN.png"),
                            line11: HandlerLineImage(image_dir + "./nitro_N-C.png"),line12: HandlerLineImage(image_dir+"./nitro_NH2+-C2.png")},
                handlelength=2, labelspacing=0.0, fontsize=36, borderpad=0.15, loc=3, 
                handletextpad=0.2, borderaxespad=0.15,bbox_to_anchor=(1.6, 0.0),prop={'size': 80})
    
    
        if element == 'H':
            line1,  = plt.plot(0,0, color='palegreen', lw=10)
            line2,  = plt.plot(0,0, color='mediumseagreen', lw=10)  
            line3, = plt.plot(0,0, color='darkseagreen', lw=10)  
            line4, = plt.plot(0,0, color='yellowgreen', lw=10)  
            line5, = plt.plot(0,0, color='darkolivegreen', lw=10)  
            line6, = plt.plot(0,0, color='pink', lw=10)  
            line7, = plt.plot(0,0, color='lightcoral', lw=10)  
            line8, = plt.plot(0,0, color='indianred', lw=10)  
            line9, = plt.plot(0,0, color='palevioletred', lw=10)  
            line10, = plt.plot(0,0, color='mediumvioletred', lw=10)  
            line11, = plt.plot(0,0, color='lightblue', lw=10)  
            line12, = plt.plot(0,0, color='lightskyblue', lw=10)  
            line13, = plt.plot(0,0, color='cornflowerblue', lw=10)  
            line14, = plt.plot(0,0, color='aquamarine', lw=10)  
            line15, = plt.plot(0,0, color='steelblue', lw=10)  
            line16, = plt.plot(0,0, color='khaki', lw=10)  
            line17, = plt.plot(0,0, color='hotpink', lw=10)  
            line18, = plt.plot(0,0, color='navajowhite', lw=10)  
            line19, = plt.plot(0,0, color='thistle', lw=10)  
            line20, = plt.plot(0,0, color='darksalmon', lw=10)  
            line21, = plt.plot(0,0, color='tomato', lw=10)  
            line22, = plt.plot(0,0, color='burlywood', lw=10) 
            
            l1 = ax4.legend([line1,line2,line3,line4,line5,
                             line6,line7,line8,line9,line10,
                             line11,line12,line13,line14,line15,
                             line16,line17,line18,line19,line20,
                             line21,line22], ["","","","","","","","","","","","","","","","","","","","","",""],
               handler_map={line1: HandlerLineImage(image_dir + "./hydro_CH3-C.png"), line2: HandlerLineImage(image_dir + "./hydro_CH2-CC.png"),
                            line3: HandlerLineImage(image_dir + "./hydro_CH-CCC.png"), line4: HandlerLineImage(image_dir + "./hydro_CH-CC.png"),
                            line5: HandlerLineImage(image_dir + "./hydro_CH-C.png"),   line6: HandlerLineImage(image_dir + "./hydro_CH3-O.png"),
                            line7: HandlerLineImage(image_dir + "./hydro_CH2-CO.png"), line8: HandlerLineImage(image_dir + "./hydro_CH-CCO.png"),
                            line9: HandlerLineImage(image_dir + "./hydro_CH-CO.png"), line10: HandlerLineImage(image_dir + "./hydro_CH-NO.png"),
                            line11: HandlerLineImage(image_dir + "./hydro_CH3-N.png"), line12: HandlerLineImage(image_dir + "./hydro_CH2-CN.png"),
                            line13: HandlerLineImage(image_dir + "./hydro_CH-CCN.png"), line14: HandlerLineImage(image_dir + "./hydro_CH-CN.png"),
                            line15: HandlerLineImage(image_dir + "./hydro_CH-NN.png"), line16: HandlerLineImage(image_dir + "./hydro_NH2-C4.png"),
                            line17: HandlerLineImage(image_dir + "./hydro_NH2-C3.png"), line18: HandlerLineImage(image_dir + "./hydro_NH-CC.png"),
                            line19: HandlerLineImage(image_dir + "./hydro_NH-C.png"),  line20: HandlerLineImage(image_dir + "./hydro_OH-C4.png"),
                            line21: HandlerLineImage(image_dir + "./hydro_OH-C3.png"), line22: HandlerLineImage(image_dir + "./hydro_OH-N.png"),
                            line1: HandlerLineImage(image_dir + "./hydro_CH3-C.png"), line2: HandlerLineImage(image_dir + "./hydro_CH2-CC.png"),
                            line3: HandlerLineImage(image_dir + "./hydro_CH-CCC.png"), line4: HandlerLineImage(image_dir + "./hydro_CH-CC.png"),
                            line5: HandlerLineImage(image_dir + "./hydro_CH-C.png"),   line6: HandlerLineImage(image_dir + "./hydro_CH3-O.png"),
                            line7: HandlerLineImage(image_dir + "./hydro_CH2-CO.png"), line8: HandlerLineImage(image_dir + "./hydro_CH-CCO.png"),
                            line9: HandlerLineImage(image_dir + "./hydro_CH-CO.png"), line10: HandlerLineImage(image_dir + "./hydro_CH-NO.png"),
                            line11: HandlerLineImage(image_dir + "./hydro_CH3-N.png"), line12: HandlerLineImage(image_dir + "./hydro_CH2-CN.png"),
                            line13: HandlerLineImage(image_dir + "./hydro_CH-CCN.png"), line14: HandlerLineImage(image_dir + "./hydro_CH-CN.png"),
                            line15: HandlerLineImage(image_dir + "./hydro_CH-NN.png"), line16: HandlerLineImage(image_dir + "./hydro_NH2-C4.png"),
                            line17: HandlerLineImage(image_dir + "./hydro_NH2-C3.png"), line18: HandlerLineImage(image_dir + "./hydro_NH-CC.png"),
                            line19: HandlerLineImage(image_dir + "./hydro_NH-C.png"),  line20: HandlerLineImage(image_dir + "./hydro_OH-C4.png"),
                            line21: HandlerLineImage(image_dir + "./hydro_OH-C3.png"), line22: HandlerLineImage(image_dir + "./hydro_OH-N.png")},
               handlelength=2, labelspacing=0.0, fontsize=36, borderpad=0.15, loc=3, 
                handletextpad=0.2, borderaxespad=0.15,bbox_to_anchor=(1.0, -1.3),prop={'size': 50})

    plt.show()

def colmark(element):
    if element == 'all':
        colors = ['lightgrey','black','blue','red']
        markers = ['.','.','+','.','.',
                   '.','.','.','.','.',
                   '.','.','*','.','.',
                   '.','.','.','*','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',]
    if element == 'H':
        markers = ['.','.','+','.','.',
                   '.','.','.','.','.',
                   '.','.','*','.','.',
                   '.','.','.','*','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',]
        
        colors=['purple','none','palegreen','pink','lightblue',
             'mediumseagreen','lightskyblue','lightcoral','none','darkseagreen',
             'cornflowerblue','indianred','none','yellowgreen','palevioletred',
             'aquamarine','none','steelblue','mediumvioletred','none',
             'darkolivegreen','none','none','khaki','hotpink',
             'navajowhite','none','none','thistle','darksalmon',
             'tomato','none','burlywood']
    if element == 'C':
        markers = ['.','.','+','.','.',
                   '.','.','.','.','.',
                   '.','.','*','.','.',
                   '.','.','.','*','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',]
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
        markers = ['.','.','+','.','.',
                   '.','.','.','.','.',
                   '.','.','*','.','.',
                   '.','.','.','*','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',]
        colors=['none','none','none','none','none',
                'none','none','none','none','none',
                'none','none','none','none','none',
                'none','none','blue','blue','blue',
                'none','none','none','none','none',
                'none','none','none','none','none',
                'none','none','none','none','none',
                'none','none','none','none','none',
                'none','none','none','none','none',
                'none','none','none','none','none',]
    if element == 'O':
        colors=['none','none','tomato','none','mediumaquamarine',
             'mediumseagreen','palegreen','none','lightblue','none',
             'lightskyblue','none','cornflowerblue','khaki','NONE',
             'none','lightpink','mistyrose','deeppink','hotpink',
             'mediumvioletred','plum','purple','NONE','mediumorchid']
        markers = ['.','.','+','.','.',
                   '.','.','.','.','.',
                   '.','.','*','.','.',
                   '.','.','.','*','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',
                   '.','.','.','.','.',]
        
    return colors, markers


def get_labels(label_file,element):
    label = pd.read_csv(label_file,delimiter=',')
    if element == 'O':
        label = {0: 'unknown',
                   1: 'HOH',
                   2: '(C)-O-(C)',
                   3: '(C)-O-(N)',
                   4: 'HO-C-(C)3',
                   5: 'HO-CH-(C)2',
                   6: 'HO-CH2-(C)',
                   7: 'HO-CH3',
                   8: 'HO-C=(C)2',
                   9: 'HOC=(O,C)',
                  10: 'HOC=(N,C)',
                  11: 'HOC=(N,O)',
                  12: 'HOC=(N,N)',
                  13: 'HO-(N)',
                  14: 'N-O-(N)',
                  15: 'O=CH2',
                  16: 'O=C-(C2)',
                  17: 'O=CH-(C)',
                  18: 'O=C-(C,N)',
                  19: 'O=C-(C,O)',
                  20: 'O=CH-(N)',
                  21: 'O=C-(N,O)',
                  22: 'O=C-(N,N)',
                  23: 'O=C-(O,O)',
                  24: 'O=CH-(O)'}
    if element == 'N':
        label = {0: 'unknown',
                 1: 'NH3',
                 2: 'NH2-C',
                 3: 'NH2-C4',
                 4: 'NH2-C3',
                 5: '',
                 6: '',
                 7: '',
                 8: '',
                 9: '',
                10: '',
                11: ''}
    if element == 'H':
        label = {0: '',
                 1: '',
                 2: '',
                 3: '',
                 4: '',
                 5: '',
                 6: '',
                 7: '',
                 8: '',
                 9: '',
                10: '',
                11: '',
                12: '',
                13: '',
                14: '',
                15: '',
                16: '',
                17: '',
                18: '',
                19: '',
                20: '',
                21: '',
                22: '',
                23: '',
                24: '',
                25: '',
                26: '',
                27: '',
                28: '',
                29: '',
                30: ''}
    return label

def extractneighborae(element,label_file):
    label = pd.read_csv(label_file,delimiter=',')
    neighbor = np.zeros((len(label),2))
    if element == 'H':
        for i in range(len(label)):
            if 1 <= int(label['Target'][i]) <= 21:
                if label['Target'][i] == 1 :
                    neighbor[i][0] = -1029.86
                    neighbor[i][1] = -13.61
                if label['Target'][i] == 2 or label['Target'][i] == 5 or label['Target'][i] == 9 or label['Target'][i] == 13 or label['Target'][i] == 20:
                    neighbor[i][0] = -1029.86
                    neighbor[i][1] = -1029.86
                if label['Target'][i] == 3 or label['Target'][i] == 7 or label['Target'][i] == 8 or label['Target'][i] == 11 or label['Target'][i] == 12 or label['Target'][i] == 14 or label['Target'][i] == 16 or label['Target'][i]==18 or label['Target'][i]==19:
                    neighbor[i][0] = -1029.86
                    neighbor[i][1] = -2042.61
                if label['Target'][i] == 4 or label['Target'][i] == 6 or label['Target'][i] == 10 or label['Target'][i] == 15 or label['Target'][i] == 17 or label['Target'][i] == 21:
                    neighbor[i][0] = -1029.86
                    neighbor[i][1] = -1485.30   
                
            if 22 <= int(label['Target'][i]) <= 28:
                if label['Target'][i] == 22:
                    neighbor[i][0] = -1485.30
                    neighbor[i][1] = -13.61
                if label['Target'][i] == 23 or label['Target'][i] == 24 or label['Target'][i] == 25 or label['Target'][i] == 28:  
                    neighbor[i][0] = -1485.30
                    neighbor[i][1] = -1029.86
                if label['Target'][i] == 26 or label['Target'][i] == 27: 
                    neighbor[i][0] = -1485.30
                    neighbor[i][1] = -1485.30   
                    
            if 29 <= int(label['Target'][i]) <= 32:
                if label['Target'][i] == 31 :
                    neighbor[i][0] = -2042.61
                    neighbor[i][1] = -13.61
                if label['Target'][i] == 29 or label['Target'][i] == 30:
                    neighbor[i][0] = -2042.61
                    neighbor[i][1] = -1029.86
                if label['Target'][i] == 32:
                    neighbor[i][0] = -2042.61
                    neighbor[i][1] = -1485.30   
                
    return neighbor