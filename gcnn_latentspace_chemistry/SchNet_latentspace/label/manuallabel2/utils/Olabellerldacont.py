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
    
    if total == 2:
        if countH == 2:
            fg_key='OHH'
            ldalabel=124
            gnucolor=16711680 
            gnumarker=13     
        if countC == 2:
            fg_key='O-CC'
            ldalabel=125
            gnucolor=16747520 
            gnumarker=7     
        if countC == 1 and countN == 1:
            fg_key='O-CN'
            ldalabel=126
            gnucolor=16776960 
            gnumarker=11   

        if countH == 1 and countC == 1:
            
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
                    fg_key='OH-C-CCC'
                    ldalabel=127
                    gnucolor=6737322 
                    gnumarker=1   
                if countC2 == 2 and countH2 == 1:
                    fg_key='OH-CH-CC'
                    ldalabel=128
                    gnucolor=3978097 
                    gnumarker=15   
                if countC2 == 1 and countH2 == 2:
                    fg_key='OH-CHH-C'
                    ldalabel = 129
                    gnucolor=10025880
                    gnumarker=9
                if countH2 == 3:
                    fg_key='OH-CHHH'
                    ldalabel = 130
                    gnucolor=14423100
                    gnumarker=5
                #label 27 & 28 is NEWWWW
#                if countO2 == 2 and countC2 ==2:
#                    label = 27
#                if countN2 == 1 and countH2 == 2:
#                    label = 28
            
            if total2 == 3:

                if countC2 == 2:
                    fg_key='OH-C-CC'
                    ldalabel = 131
                    gnucolor=11393254
                    gnumarker=2
                if countC2 == 1 and countO2 == 2:
                    fg_key='OH-C-CO'
                    ldalabel = 132
                    gnucolor=6266528
                    gnumarker=5
                if countC2 == 1 and countN2 == 1:
                    fg_key='OH-C-CN'
                    ldalabel = 133
                    gnucolor=8900346
                    gnumarker=13
                if countO2 == 2 and countN2 == 1:
                    fg_key='OH-C-NO'
                    ldalabel = 134
                    gnucolor=4286945
                    gnumarker=13
                if countN2 == 2:
                    fg_key='OH-C-NN'
                    ldalabel = 135
                    gnucolor=6591981
                    gnumarker=3
                #label 26 is NEW
#                if countC2 == 1 and countH2 == 1:
#                    label = 26
                #label 29 is NEW
#                if countO2 == 2 and countH2 == 1:
#                    label = 29
            
                

                
        if countH == 1 and countN == 1:
            fg_key='OH-N'
            ldalabel = 136
            gnucolor=15787660
            gnumarker=7
        if countN == 2:
            fg_key='O-NN'
            ldalabel = 137
            gnucolor=16766720
            gnumarker=5
    if total == 1:


        if countC == 1:
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

            #for fgtransform
            if total2 == 4:
                if countH2 ==3 and countC2 == 1:
                    ldalabel = 138
                    gnucolor = 2142890
                    gnumarker = 15
            if total2 == 3: 
                if countH2 == 2:
                    fg_key='O-C-HH'
                    ldalabel = 139
                    gnucolor=11674146
                    gnumarker=15
                if countC2 == 2:
                    fg_key='O-C-CC'
                    ldalabel = 140
                    gnucolor=16758465
                    gnumarker=11
                if countC2 == 1 and countH2 == 1:
                    fg_key='O-CH-C'
                    ldalabel = 141
                    gnucolor=16770273
                    gnumarker=1
                if countC2 == 1 and countN2 == 1:
                    fg_key='O-C-CN'
                    ldalabel = 142
                    gnucolor=16716947
                    gnumarker=3
                if countC2 == 1 and countO2 == 2:
                    fg_key='O-C-CO'
                    ldalabel = 143
                    gnucolor=16738740
                    gnumarker=15
                if countN2 == 1 and countH2 == 1:
                    fg_key='O-CH-N'
                    ldalabel=144
                    gnucolor=13047173
                    gnumarker=3
                if countN2 == 1 and countO2 == 2:
                    fg_key='O-C-NO'
                    ldalabel=145
                    gnucolor=14524637
                    gnumarker=5
                #label = 25 IS NEW
#                if countN2 == 2 and countO2 == 1:
#                    label = 25
                if countN2 == 2:
                    fg_key='O-C-NN'
                    ldalabel=146
                    gnucolor=8388736
                    gnumarker=5
                if countO2 == 3:
                    fg_key='O-C-OO'
                    ldalabel=147
                    gnucolor=14381203
                    gnumarker=9
                if countO2 == 2 and countH2 == 1:
                    fg_key='O-CH-O'
                    ldalabel=148
                    gnucolor=12211667
                    gnumarker=7
        if countN == 1:
            fg_key='O-N'
            ldalabel=149
            gnucolor=10824234
            gnumarker=7

    if ldalabel == 999:
        print('999')
        error_message = 'this O-type functional group is unknown, molecule_id = %s atom_index =  %s, H, %s %s, C, %s %s, N, %s %s, O, %s %s, F, %s %s' %(each_molecule,each_atom,countH,countH2,countC,countC2,countN,countN2,countO,countO2,countF,countF2)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker