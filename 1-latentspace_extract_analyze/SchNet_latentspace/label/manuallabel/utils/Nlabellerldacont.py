from label.manuallabel2.utils import utils

def label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom):
    countC = 0
    countH = 0
    countO = 0
    countN = 0
    countF = 0
    for each_neighbor in range(len(connections_to_atom)):

        for line_index, each_line in enumerate(xyz_file_read.split('\n')):
            if line_index == int(connections_to_atom[each_neighbor])+1:
                if '1' in each_line[0:1]:
                    countH = countH+1
                if '6' in each_line[0:1]:
                    countC = countC+1
                if '7' in each_line[0:1]:
                    countN = countN+1
                if '8' in each_line[0:1]:
                    countO = countO+1
                if '9'in each_line[0:1]:
                    countF = countF+1
                    
    total = countC + countH + countO + countN + countF
    ldalabel = 999

    if total == 3:
        if countH == 3:
            fg_key='NHHH'
            ldalabel=96
            gnucolor=16711680 
            gnumarker=13     
        if countH == 2 and countC == 1:


            #locate the carbon neighbor
            check_if_hydrogen = True
            while check_if_hydrogen == True:
                check_neighbor = min(connections_to_atom)
                element_present = utils.check_element(int(check_neighbor),xyz_file_read)
                if element_present == 'H':
                    check_if_hydrogen = True
                    index = connections_to_atom.index(int(check_neighbor))
                    connections_to_atom.pop(index)


                else:
                    carbon_neighbor = check_neighbor
                    check_if_hydrogen = False
            
            #get the neighbors of this carbon
            carbon_neighbor_neighbors = utils.find_connections(number_atoms,carbon_neighbor,mol_file_read)
            total2,countC2,countH2,countO2,countN2, countF2 = utils.count_nn(carbon_neighbor_neighbors,xyz_file_read)

            if total2 == 4: 
                if countC2 == 3:
                    fg_key='NHH-C-CCC'
                    ldalabel=97
                    gnucolor=16753920 
                    gnumarker=7  
                if countH2 == 1 and countC2 == 2:
                    fg_key='NHH-CH-CC'
                    ldalabel=98
                    gnucolor=16776960 
                    gnumarker=11 
                if countH2 == 2 and countC2 == 1:
                    fg_key='NHH-CHH-C'
                    ldalabel=99
                    gnucolor=6737322 
                    gnumarker=1 

            if total2 == 3:
                if countC2 == 2:
                    fg_key='NHH-C-CC'
                    ldalabel=100
                    gnucolor=3978097 
                    gnumarker=15
                if countC2 == 1 and countO2 == 1:
                    fg_key='NHH-C-CO'
                    ldalabel=101
                    gnucolor=10025880 
                    gnumarker=9
                if countO2 == 2:
                    fg_key='NHH-C-OO'
                    ldalabel=102
                    gnucolor=10145074
                    gnumarker=5
                if countO2 == 1 and countN2 == 2:
                    fg_key='NHH-C-NO'
                    ldalabel=103
                    gnucolor=11393254
                    gnumarker=15
                if countN2 == 3: 
                    fg_key='NHH-C-NN'
                    ldalabel=104
                    gnucolor=6266528
                    gnumarker=3
                if countN2 == 2 and countC2 == 1:
                    fg_key='NHH-C-CN'
                    ldalabel= 105
                    gnucolor=8900346
                    gnumarker=13
                if countN2 == 2 and countH2 == 1:
                    fg_key='NHH-CH-N'
                    ldalabel= 106
                    gnucolor=4286945
                    gnumarker=13
                if countO2 == 1 and countH2 ==1:
                    fg_key='NHH-CH-O'
                    ldalabel= 107
                    gnucolor=6591981
                    gnumarker=3
    
        if countH == 1 and countC == 2:
            fg_key='NH-CC'
            ldalabel= 108
            gnucolor=15787660
            gnumarker=7
        if countH == 1 and countN == 2:
            fg_key='NH-NN'
            ldalabel= 109
            gnucolor=16766720
            gnumarker=5
        if countH ==1 and countN == 1 and countC ==1:
            fg_key='NH-CN'
            ldalabel= 110
            gnucolor=13468991
            gnumarker=15
        if countC == 3:
            fg_key='N-CCC'
            ldalabel= 111
            gnucolor=16758465
            gnumarker=11
        if countC == 2 and countN == 1:
            fg_key='N-CCN'
            ldalabel= 112
            gnucolor=16770273
            gnumarker=1
        if countC == 1 and countN ==2:
            fg_key='N-CNN'
            ldalabel= 113
            gnucolor=16716947
            gnumarker=3
    if total == 2:
        if countC == 2:
            fg_key='N-CC'
            ldalabel= 114
            gnucolor=16738740
            gnumarker=15
        if countC == 1 and countH == 1:
            fg_key='NH-C'
            ldalabel= 115
            gnucolor=13047173
            gnumarker=1
        if countC == 1 and countO == 1:
            fg_key='N-CO'
            ldalabel= 116
            gnucolor=14524637
            gnumarker=13
        if countC == 1 and countN == 1:
            fg_key='N-CN'
            ldalabel= 117
            gnucolor=8388736
            gnumarker=5
        if countO == 1 and countN == 1:
            fg_key='N-NO'
            ldalabel= 118
            gnucolor=14381203
            gnumarker=9 
        if countN == 2:
            fg_key='N-NN'
            ldalabel= 119
            gnucolor=12211667
            gnumarker=7 
    if total == 1:
        if countC == 1:
            fg_key='N-C'
            ldalabel= 120
            gnucolor=13808780
            gnumarker=5 
    if total == 4:
        if countH == 3 and countC == 1:
            fg_key='NHHH-C'
            ldalabel= 121
            gnucolor=13458524
            gnumarker=7 
        if countH == 2 and countC == 2:
            fg_key='NHH-CC'
            ldalabel= 122
            gnucolor=16752762
            gnumarker=7
        if countC == 3 and countH == 1:
            fg_key='NH-CCC'
            ldalabel= 123
            gnucolor=15761536
            gnumarker=9
#        if countC == 4:
#            label = 32 #NEW

    if ldalabel == 999:
        print('999')
        error_message = 'this N-type functional group is unknown, molecule_id = %s, atom_index = %s, H, %s, C, %s, N, %s, O, %s, F, %s' %(each_molecule,each_atom,countH,countC,countN,countO,countF)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker