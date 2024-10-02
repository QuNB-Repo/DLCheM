import matplotlib
import matplotlib.pyplot as plt

import numpy as np
from numpy import genfromtxt, savetxt


def plot2d(data_filepath,label_filepath,save_filepath,labelid,component1,component2,color,colors,marker,markers,x_label,y_label):

    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')
    print(data)
    fig,ax = plt.subplots(figsize=(10,6),dpi=500,nrows=1,ncols=1)

    if marker == True:  
        for i in range(len(data)):
            label = genfromtxt(label_filepath,delimiter=' ',encoding='utf-8-sig')
            j = label[i][0]
            j = int(j)
            ax.scatter(data[i][component1],data[i][component2],color=colors[j],marker=markers[j])
    elif color == True:
        label = genfromtxt(label_filepath,delimiter=' ',encoding='utf-8-sig')
        ax.scatter(data[:,component1],data[:,component2],c=label[:,0],cmap=matplotlib.colors.ListedColormap(colors))
    else: 
        ax.scatter(data[:,component1],data[:,component2])

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(save_filepath)
    plt.show()

    for i in range(1,len(labelid)):
        print(i, labelid[i])

    print(save_filepath)

def lines(data_filepath,n_data,divide,save_filepath,x_label,y_label):

    data = genfromtxt(data_filepath,delimiter=',',encoding='utf-8-sig')

    fig, ax = plt.subplots(figsize=(10,6),dpi=100,nrows=1,ncols=1)

    #for each pca data point in the first 138, 
    for i in range(int(n_data)):
        x1 = data[i][0]
        y1 = data[i][1]

        x2 = data[i+int(divide)][0]
        y2 = data[i+int(divide)][1]

        coordinate1 = [x1,y1]
        coordinate2 = [x2,y2]

        coordinate = np.array([coordinate1,coordinate2])
        print(coordinate)

        ax.plot(coordinate[:,0],coordinate[:,1])

        plt.show()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.savefig(save_filepath)


    # coordinates of each 

def gnuplot(data_file,label_file,save_file):
    data = genfromtxt(data_file,delimiter=',')
    label = genfromtxt(label_file,delimiter=',')
    data = np.column_stack((data,label[:,0]))
    savetxt(save_file,data)

def exp2d(xrange,yrange,data_file,label_file,save_file):
    data = genfromtxt(data_file,delimiter=',',encoding='utf-8-sig')
    label = genfromtxt(label_file,delimiter=',',encoding='utf-8-sig')

    new_data = []
    for i in range(len(data)):
        if  xrange[0] < data[i][0] < xrange[1]:
            if yrange[0] < data[i][1] < yrange[1]:
                new_data.append([data[i][0],data[i][1],label[i][0],label[i][1]])
                

    savetxt(save_file,new_data,delimiter=',')

def colmark(element):
    if element == 'all':
        colors = ['pink','silver','lightskyblue','tomato','violet']
        markers = ['P','D','s','o','*']
    if element == 'H':
        markers=['D','o','s','P','^',
                 'P','*','X','o','v',
                 'X','*','p','v','s',
                 'p','P','*','D','s',
                 'o','X','s','*','s',
                 '^','P','*','o','v',
                 '>','X','D','P','p',
                 '*','<','>','o','<',
                 's','*','X','D','o',
                 'X','P','v','D','^',
                 '*']

        colors=['red','palegreen','yellowgreen','limegreen','paleturquoise',
                'mediumseagreen','darkmagenta','steelblue','darkslategrey','darkcyan',
                'darkslateblue','midnightblue','darkviolet','olive','royalblue',
                'orchid','blueviolet','cadetblue','firebrick','green',
                'thistle','crimson','palegoldenrod','peru','navajowhite',
                'yellow','gold','goldenrod','darkkhaki','tomato',
                'darkgoldenrod','sandybrown','orange','darkorange','chocolate',
                'khaki','purple','brown','grey','peachpuff',
                'lightsalmon','coral','salmon','lightpink','hotpink',
                'deeppink','palevioletred','magenta','violet','maroon',
                'moccasin']
        
        label = {0: 'CH4, ' + colors[0] + ', ' + markers[0],
                 1: 'CH3-C, '+ colors[1] + ', ' + markers[1],
                 2: 'CH3-O, ' + colors[2] + ', ' + markers[2],
                 3: 'CH3-N, ' + colors[3] + ', ' + markers[3],
                 4: 'C-CH2-C, ' + colors[4] + ', ' + markers[4],
                 5: 'N-CH2-C, ' + colors[5] + ', ' + markers[5],
                 6: 'O-CH2-C, ' + colors[6] + ', ' + markers[6],
                 7: 'O-CH2-O, '  + colors[7] + ', ' + markers[7],
                 8: 'C-CH-C(2), ' +  colors[8] + ', ' + markers[8],
                 9: 'N-CH-C(2), ' +  colors[9] + ', ' + markers[9],
                10: 'O-CH-C(2), ' + colors[10] + ', ' + markers[10],
                11: 'O-CH-CO, ' + colors[11] + ', ' + markers[11],
                12: 'C-CH-C, ' + colors[12] + ', ' + markers[12],
                13: 'C-CH-O, ' + colors[13] + ', ' + markers[13],
                14: 'C-CH-N, ' + colors[14] + ', ' + markers[14],
                15: 'O-CH-O, '+ colors[15] + ', ' + markers[15],
                16: 'N-CH-N, '+ colors[16] + ', ' + markers[16],
                17: 'N-CH-O, '+ colors[17] + ', ' + markers[17] ,
                18: 'H2C=O, ' + colors[18] + ', ' + markers[18],
                19: 'HC=-C, '+ colors[19] + ', ' + markers[19],
                20: 'HC=-N, '+ colors[20] + ', ' + markers[20],
                21: 'NH3, '+ colors[21] + ', ' + markers[21],
                22: 'NH2-C-C(3), ' + colors[22] + ', ' + markers[22],
                23: 'NH2-CH-C(2), ' + colors[23] + ', ' + markers[23],
                24: 'NH2-CH2-C, ' + colors[24] + ', ' + markers[24],
                25: 'NH2-C=CC, '+ colors[25] + ', ' + markers[25],
                26: 'NH2-C=CN, '+ colors[26] + ', ' + markers[26],
                27: 'NH2-C=CO, '+ colors[27] + ', ' + markers[27],
                28: 'NH2-C=O(2), ' + colors[28] + ', ' + markers[28],
                29: 'NH2-CH=O, ' + colors[29] + ', ' + markers[29],
                30: 'NH2-C=NO, ' + colors[30] + ', ' + markers[30],
                31: 'NH2-C=N(2), ' + colors[31] + ', ' + markers[31],
                32: 'C-NH-C, ' + colors[32] + ', ' + markers[32],
                33: 'N-NH-C, '+ colors[33] + ', ' + markers[33],
                34: 'N-NH-N, '+ colors[34] + ', ' + markers[34],
                35: 'HN=-C, '+ colors[35] + ', ' + markers[35],
                36: 'NH3-C-' + colors[36] + ',' + markers[36],
                37: 'C-NH2-C' + colors[37] + ',' + markers[37],
                38: 'OH-CH3' + colors[38] + ',' + markers[38], 
                39: 'OH-C-C(3), ' + colors[39] + ', ' + markers[39],
                40: 'OH-CH-C(2), ' + colors[40] + ', ' + markers[40],
                41: 'OH-C-C(2)N, ' + colors[41] + ', ' + markers[41],
                42: 'OH-C-C(2)O, ' + colors[42] + ', ' + markers[42],
                43: 'OH-CH2-C, ' + colors[43] + ', ' + markers[43],
                44: 'OH-C=CC, '+ colors[44] + ', ' + markers[44],
                45: 'OH-C=CN, '+ colors[45] + ', ' + markers[45],
                46: 'OH-C=CO, '+ colors[46] + ', ' + markers[46],
                47: 'OH-C=NO, ' + colors[47] + ', ' + markers[47],
                48: 'OH-C=N(2), ' + colors[48] + ', ' + markers[48],
                49: 'H-O-H, '+ colors[49] + ', ' + markers[49],
                50: 'H-O-N, '+ colors[50] + ', ' + markers[50],
                }
        
    if element == 'C':
         markers=['D','o','s','P','^',
                 'P','*','X','o','v',
                 'X','*','>','v','s',
                 'p','P','*','D','s',
                 'o','D','s','*','s',
                 '^','P','*','o','v',
                 '>','X','D','P','p',
                 '*','<','>','o','<',
                 's','*']

         colors=['forestgreen','palegreen','yellowgreen','limegreen','red',
                'paleturquoise','darkmagenta','steelblue','purple','darkcyan',
                'cornflowerblue','midnightblue','skyblue','darkslateblue','royalblue',
                'firebrick','mediumvioletred','lightseagreen','darkseagreen','green',
                'crimson','dodgerblue','darkslategrey','mediumaquamarine','lightsteelblue',
                'yellow','gold','magenta','hotpink','deeppink',
                'goldenrod','mediumpurple','orange','mediumslateblue','palevioletred',
                'pink','rosybrown','indianred','olive','orangered',
                'burlywood','peru']
        
         label = {0: 'C-(C4), ' + colors[0] + ', ' + markers[0],
                 1: 'CH-(C3), '+ colors[1] + ', ' + markers[1],
                 2: 'CH2-(C2), ' + colors[2] + ', ' + markers[2],
                 3: 'CH3-C, ' + colors[3] + ', ' + markers[3],
                 4: 'CH4, ' + colors[4] + ', ' + markers[4],
                 5: 'CH2-(CO), ' + colors[5] + ', ' + markers[5],
                 6: 'CH2-(CN), ' + colors[6] + ', ' + markers[6],
                 7: 'CH3-(O), ' +  colors[7] + ', ' + markers[7],
                 8: 'CH3-(N), ' +  colors[8] + ', ' + markers[8],
                 9: 'CH-(C2O), ' + colors[9] + ', ' + markers[9],
                10: 'CH-(CO2), ' + colors[10] + ', ' + markers[10],
                11: 'CH-(C2N), ' + colors[11] + ', ' + markers[11],
                12: 'C-(C3O), ' + colors[12] + ', ' + markers[12],
                13: 'C-(C3N), '+ colors[13] + ', ' + markers[13],
                14: 'CH2-(O2), '+ colors[14] + ', ' + markers[14],
                15: 'CF4, '+ colors[15] + ', ' + markers[15],
                16: 'C-(CF3), '+ colors[16] + ', ' + markers[16],
                17: 'C-(C2O2), '+ colors[17] + ', ' + markers[17],
                18: 'C-(C3), ' + colors[18] + ', ' + markers[18],
                19: 'CH(C2), ' + colors[19] + ', ' + markers[19],
                20: 'CH2-O, '+ colors[20] + ', ' + markers[20],
                21: 'CH-(CO), '+ colors[21] + ', ' + markers[21],
                22: 'C-(C2O), '+ colors[22] + ', ' + markers[22],
                23: 'C-(CO2), ' + colors[23] + ', ' + markers[23],
                24: 'CH-(O2), ' + colors[24] + ', ' + markers[24],
                25: 'CH-(NO), ' + colors[25] + ', ' + markers[25],
                26: 'C-(CNO), ' + colors[26] + ', ' + markers[26],
                27: 'C-(CN2), ' + colors[27] + ', ' + markers[27],
                28: 'CH-(CN), '+ colors[28] + ', ' + markers[28],
                29: 'C-(N2O), '+ colors[29] + ', ' + markers[29],
                30: 'C-(NO2), '+ colors[30] + ', ' + markers[30],
                31: 'C-(N3)-' + colors[31] + ',' + markers[31],
                32: 'C-(O3)' + colors[32] + ',' + markers[32],
                33: 'C-(C2N)' + colors[33] + ',' + markers[33], 
                34: 'CH-(N2), ' + colors[34] + ', ' + markers[34],
                35: 'C-(C2F), ' + colors[35] + ', ' + markers[35],
                36: 'C-(CNF), ' + colors[36] + ', ' + markers[36],
                37: 'C-(N2F), ' + colors[37] + ', ' + markers[37],
                38: 'C-(C2), ' + colors[38] + ', ' + markers[38],
                39: 'CH-C, '+ colors[39] + ', ' + markers[39],
                40: 'CH-N, '+ colors[40] + ', ' + markers[40],
                41: 'C-(CN), '+ colors[41] + ', ' + markers[41],
                }

    if element == 'N':
        colors=['red','orange','yellow','mediumaquamarine','mediumseagreen',
                'palegreen','yellowgreen','lightblue','cadetblue','lightskyblue',
                'royalblue','cornflowerblue','khaki','gold','peru',
                'lightpink','mistyrose','deeppink','hotpink','mediumvioletred',
                'plum','purple','palevioletred','mediumorchid','tan',
                'indianred','lightsalmon','lightcoral']
        markers = ['D','o','v','P','p',
                   '^','s','X','*','<',
                   'D','*','o','s','>',
                   'v','P','*','p','P',
                   'D','s','^','o','s',
                   'o','^','X']
        
        label = {0: 'NH3, '+ colors[0] +', ' + markers[0], 
                 1: 'NH2-C-(C3), ' + colors[1] + ', ' + markers[1],
                 2: 'NH2-CH-(C2), ' + colors[2] + ', ' + markers[2],
                 3: 'NH2-CH2-C, ' + colors[3] + ', ' + markers[2],
                 4: 'NH2-C-(C2), ' + colors[4] + ', ' + markers[4],
                 5: 'NH2-C-(CO), ' + colors[5] + ', ' + markers[5],
                 6: 'NH2-C-(OO), ' + colors[6] + ', ' + markers[6],
                 7: 'NH2-C-(NO), ' + colors[7] + ', ' + markers[7],
                 8: 'NH2-C-(NN), ' + colors[8] + ', ' + markers[8],
                 9: 'NH2-C-(CN), ' + colors[9] + ', ' + markers[9],
                 10: 'NH2-CH-N, ' + colors[10] + ', ' + markers[10] ,
                 11: 'NH2-CH-O, ' + colors[11] + ', ' + markers[11],
                 12: 'C-NH-C, ' + colors[12] + ', ' + markers[12],
                 13: 'N-NH-N, '+ colors[13] + ', ' + markers[13],
                 14: 'N-NH-C, ' + colors[14] + ', ' + markers[14],
                 15: '(C)2-N-C, ' + colors[15] + ', ' + markers[15],
                 16: '(C)2-N-N, ' + colors[16] + ', ' + markers[16],
                 17: 'C-N-(N)2, ' + colors[17] + ', ' + markers[17],
                 18: 'C-N-C, ' + colors[18] + ', '  + markers[18],
                 19: 'NH-C, ' + colors[19] + ', ' + markers[19],
                 20: 'C-N-O, ' + colors[20] + ', ' + markers[20],
                 21: 'C-N-N, ' + colors[21] + ', ' + markers[21],
                 22: 'N-N-O, ' + colors[22] + ', ' + markers[22],
                 23: 'N-N-N, ' + colors[23] + ', ' + markers[23],
                 24: 'N-C, ' + colors[24] + ', ' + markers[24],
                 25: 'NH3-C, ' + colors[25] + ', ' + markers[25],
                 26: 'C-NH2-C, ' + colors[26] + ', ' + markers[26],
                 27: '(C)3-NH, ' + colors[27] + ', ' + markers[27],
                 }
        
    if element == 'O':
        colors=['red','orange','yellow','mediumaquamarine','mediumseagreen',
                'palegreen','crimson','lightblue','cadetblue','lightskyblue',
                'royalblue','cornflowerblue','khaki','gold','firebrick',
                'lightpink','mistyrose','deeppink','hotpink','mediumvioletred',
                'plum','purple','palevioletred','mediumorchid','brown']
        markers = ['D','o','v','P','p',
                   '^','s','X','*','<',
                   'D','*','o','s','>',
                   'v','P','*','p','*',
                   '*','s','^','o','o',
                   ]
        
        label = {0: 'HOH, '+ colors[0] +', ' + markers[0],
                 1: '(C)-O-(C), '+ colors[1] +', ' + markers[1],
                 2: '(C)-O-(N), '+ colors[2] +', ' + markers[2],
                 3: 'HO-C-(C)3-ringed+non-ringed, '+ colors[3] +', ' + markers[3],
                 4: 'HO-CH-(C)2, '+ colors[4] +', ' + markers[4],
                 5: 'HO-CH2-(C), '+ colors[5] +', ' + markers[5],
                 6: 'HO-CH3, '+ colors[6] +', ' + markers[6],
                 7: 'HO-C=(C)2, '+ colors[7] +', ' + markers[7],
                 8: 'HOC=(O,C*), '+ colors[8] +', ' + markers[8],
                 9: 'HOC=(N,C*), '+ colors[9] +', ' + markers[9],
                10: 'HOC=(N*,O), '+ colors[10] +', ' + markers[10],
                11: 'HOC=(N,N), '+ colors[11] +', ' + markers[11],
                12: 'HO-(N), '+ colors[12] +', ' + markers[12],
                13: 'N-O-(N), '+ colors[13] +', ' + markers[13],
                14: 'O=CH2, '+ colors[14] +', ' + markers[14],
                15: 'O=C-(C2), '+ colors[15] +', ' + markers[15],
                16: 'O=CH-(C), '+ colors[16] +', ' + markers[16],
                17: 'O=C-(C,N)-ringed+non-ringed, '+ colors[17] +', ' + markers[17],
                18: 'O=C-(C,O)-ringed+non-ringed, '+ colors[18] +', ' + markers[18],
                19: 'O=CH-(N), '+ colors[19] +', ' + markers[19],
                20: 'O=C-(N,O)-ringed+non-ringed, '+ colors[20] +', ' + markers[20],
                21: 'O=C-(N,N)-ringed+non-ringed, '+ colors[21] +', ' + markers[21],
                22: 'O=C-(O,O), '+ colors[22] +', ' + markers[22],
                23: 'O=CH-(O), '+ colors[23] +', ' + markers[23],
                24: 'O=N-, '+ colors[24] +', ' + markers[24],
                  }        

        
    return colors, markers, label
