from label.manuallabel2.utils import utils

def label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom):
    countC = 0
    countH = 0
    countO = 0
    countN = 0
    countF = 0
    countP = 0
    countS = 0
    countCl = 0
    countBr = 0
    countI = 0
    for each_neighbor in range(len(connections_to_atom)):

        for line_index, each_line in enumerate(xyz_file_read.split('\n')):
            if line_index == int(connections_to_atom[each_neighbor])+1:
                if '1 ' in each_line[0:2]:
                    countH = countH+1
                if '6 ' in each_line[0:2]:
                    countC = countC+1
                if '7 ' in each_line[0:2]:
                    countN = countN+1
                if '8 ' in each_line[0:2]:
                    countO = countO+1
                if '9 'in each_line[0:2]:
                    countF = countF+1
                if '15' in each_line[0:2]:
                    countP = countP+1
                if '16' in each_line[0:2]:
                    countS = countS+1
                if '17' in each_line[0:2]:
                    countCl = countCl+1
                if '35' in each_line[0:2]:
                    countBr = countBr + 1
                if '53' in each_line[0:2]:
                    countI = countI + 1
                    
    total = countC + countH + countO + countN + countF + countP + countS + countCl + countBr + countI
    ldalabel = 999
    
    if total == 2:
        if countH == 2:
            fg_key='OHH'
            ldalabel=0
            gnucolor=16711680 
            gnumarker=13 
            hexcolor='#FF0000'    
        if countC == 2:
            fg_key='O-CC'
            ldalabel=1
            gnucolor=16747520 
            gnumarker=7 
            hexcolor='#FF8C00'      
        if countC == 1 and countN == 1:
            fg_key='O-CN'
            ldalabel=2
            gnucolor=16776960 
            gnumarker=11 
            hexcolor='#FFFF00'    

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
                    ldalabel=3
                    gnucolor=6737322 
                    gnumarker=1 
                    hexcolor='#66CDAA'    
                if countC2 == 2 and countH2 == 1:
                    fg_key='OH-CH-CC'
                    ldalabel=4
                    gnucolor=3978097 
                    gnumarker=15 
                    hexcolor='#3CB371'    
                if countC2 == 1 and countH2 == 2:
                    fg_key='OH-CHH-C'
                    ldalabel = 5
                    gnucolor=10025880
                    gnumarker=9 
                    hexcolor='#98FB98' 
                if countH2 == 3:
                    fg_key='OH-CHHH'
                    ldalabel = 6
                    gnucolor=14423100
                    gnumarker=5 
                    hexcolor='#DC143C' 

            if total2 == 3:

                if countC2 == 2:
                    fg_key='OH-C-CC'
                    ldalabel = 7
                    gnucolor=11393254
                    gnumarker=2 
                    hexcolor='#ADD8E6' 
                if countC2 == 1 and countO2 == 2:
                    fg_key='OH-C-CO'
                    ldalabel = 8
                    gnucolor=6266528
                    gnumarker=5 
                    hexcolor='#5F9EA0' 
                if countC2 == 1 and countN2 == 1:
                    fg_key='OH-C-CN'
                    ldalabel = 9
                    gnucolor=8900346
                    gnumarker=13 
                    hexcolor='#87CEFA' 
                if countO2 == 2 and countN2 == 1:
                    fg_key='OH-C-NO'
                    ldalabel = 10
                    gnucolor=4286945
                    gnumarker=13 
                    hexcolor='#4169E1' 
                if countN2 == 2:
                    fg_key='OH-C-NN'
                    ldalabel = 11
                    gnucolor=6591981
                    gnumarker=3 
                    hexcolor='#6495ED' 

                
        if countH == 1 and countN == 1:
            fg_key='OH-N'
            ldalabel = 12
            gnucolor=15787660
            gnumarker=7 
            hexcolor='#F0E68C' 
        if countN == 2:
            fg_key='O-NN'
            ldalabel = 13
            gnucolor=16766720
            gnumarker=5
            hexcolor='#FFD700' 
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


            if total2 == 3: 
                if countH2 == 2:
                    fg_key='O-C-HH'
                    ldalabel = 14
                    gnucolor=11674146
                    gnumarker=15
                    hexcolor='#B22222' 
                if countC2 == 2:
                    fg_key='O-C-CC'
                    ldalabel = 15
                    gnucolor=16758465
                    gnumarker=11
                    hexcolor='#FFB6C1' 
                if countC2 == 1 and countH2 == 1:
                    fg_key='O-CH-C'
                    ldalabel = 16
                    gnucolor=16770273
                    gnumarker=1
                    hexcolor='#FFE4E1' 
                if countC2 == 1 and countN2 == 1:
                    fg_key='O-C-CN'
                    ldalabel = 17
                    gnucolor=16716947
                    gnumarker=3
                    hexcolor='#FF1493' 
                if countC2 == 1 and countO2 == 2:
                    fg_key='O-C-CO'
                    ldalabel = 18
                    gnucolor=16738740
                    gnumarker=15
                    hexcolor='#FF69B4' 
                if countN2 == 1 and countH2 == 1:
                    fg_key='O-CH-N'
                    ldalabel=19
                    gnucolor=13047173
                    gnumarker=3
                    hexcolor='#C71585' 
                if countN2 == 1 and countO2 == 2:
                    fg_key='O-C-NO'
                    ldalabel=20
                    gnucolor=14524637
                    gnumarker=5
                    hexcolor='#DDA0DD' 
                if countN2 == 2:
                    fg_key='O-C-NN'
                    ldalabel=21
                    gnucolor=8388736
                    gnumarker=5
                    hexcolor='#800080' 
                if countO2 == 3:
                    fg_key='O-C-OO'
                    ldalabel=22
                    gnucolor=14381203
                    gnumarker=9
                    hexcolor='#DB7093' 
                if countO2 == 2 and countH2 == 1:
                    fg_key='O-CH-O'
                    ldalabel=23
                    gnucolor=12211667
                    gnumarker=7
                    hexcolor='#BA55D3' 
        if countN == 1:
            fg_key='O-N'
            ldalabel=24
            gnucolor=10824234
            gnumarker=7
            hexcolor='#A52A2A' 
    
    #NEW because of IGC50!!!!!
    if total == 2:
        if countS == 1 and countC == 1:
            fg_key = 'S-O-C'
            ldalabel=25
            gnucolor=12092939
            gnumarker=3
            hexcolor='#B8860B' 
        if countP == 1 and countC == 1:
            fg_key = 'P-O-C'
            ldalabel=26
            gnucolor=13468991
            gnumarker=1
            hexcolor='#CD853F' 

    if total == 1:
        if countS == 1:
            fg_key = 'O=S'
            ldalabel=27
            gnucolor=8421376
            gnumarker=7
            hexcolor='#808000' 
        if countP == 1:
            fg_key = 'O=P'
            ldalabel=28
            gnucolor=13808780
            gnumarker=5
            hexcolor='#D2B48C' 




    if ldalabel == 999:
        print('999')
        error_message = 'this O-type functional group is unknown, molecule_id = %s atom_index =  %s, H, %s %s, C, %s %s, N, %s %s, O, %s %s, F, %s %s' %(each_molecule,each_atom,countH,countH2,countC,countC2,countN,countN2,countO,countO2,countF,countF2)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker,hexcolor