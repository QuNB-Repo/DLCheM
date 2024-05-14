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
            if total2 == 4:
                if countC2 == 3:
                    fg_key='F-C-CCC'
                    ldalabel=0
                    gnucolor=2263842 
                    gnumarker=13    
                if countH2 == 1 and countC2 == 2:
                    fg_key='F-CH-CC'
                    ldalabel=1
                    gnucolor=10025880
                    gnumarker=7
                if countH2 == 2 and countC2 == 1:
                    fg_key='F-CHH-C'
                    ldalabel=2
                    gnucolor=10145074
                    gnumarker=5
                if countH2 == 4: 
                    fg_key='F-CHHH'
                    ldalabel=3
                    gnucolor=16711680
                    gnumarker=9
                if countH2 == 2 and countO2 == 1:
                    fg_key='F-CHH-O'
                    ldalabel=4
                    gnucolor=11529966
                    gnumarker=1
                if countH2 == 2 and countN2 == 1:
                    fg_key='F-CHH-N'
                    ldalabel=7
                    gnucolor=9109643
                    gnumarker=3
                if countH2 == 2 and countO2 == 1:
                    fg_key='F-CHH-O'
                    ldalabel=8
                    gnucolor=4620980
                    gnumarker=2
                if countH2 == 2 and countN2 == 1:
                    fg_key='F-CHH-N'
                    ldalabel=9
                    gnucolor=8388736
                    gnumarker=7
                if countH2 == 1 and countC2 == 1 and countO2 == 1:
                    fg_key='F-CH-CO'
                    ldalabel=10
                    gnucolor=35723
                    gnumarker=11
                if countH2 == 1 and countO2 == 2:
                    fg_key='F-CH-OO'
                    ldalabel=11
                    gnucolor=6591981
                    gnumarker=2
                if countH2 == 1 and countC2 == 1 and countN2 == 1:
                    fg_key='F-CH-CN'
                    ldalabel=12
                    gnucolor=1644912
                    gnumarker=3
                if countC2 == 2 and countO2 == 1:
                    fg_key='F-C-CCO'
                    ldalabel=13
                    gnucolor=8900331
                    gnumarker=15
                if countC2 == 2 and countN2 == 1:
                    fg_key='F-C-CCN'
                    ldalabel=14
                    gnucolor=4734347
                    gnumarker=11
                if countH2 == 1 and countO2 == 2:
                    fg_key='F-CH-OO'
                    ldalabel=15
                    gnucolor=4286945
                    gnumarker=5
                if countF2 == 4: 
                    fg_key='CF4'
                    ldalabel=16
                    gnucolor=11674146
                    gnumarker=15
                if countF2 == 3 and countC2 == 1:
                    fg_key='F-C-CFF'
                    ldalabel=17
                    gnucolor=13047173
                    gnumarker=1
                if countC2 == 1 and countO2 == 2:
                    fg_key='F-C-COO'
                    ldalabel=18
                    gnucolor=2142890
                    gnumarker=3
            if total2 == 3:
                if countC2 == 2:
                    fg_key='F-C-CC'
                    ldalabel=19
                    gnucolor=9419919
                    gnumarker=13
                if countH2 == 1 and countC2 == 1:
                    fg_key='F-CH-C'
                    ldalabel=20
                    gnucolor=32768
                    gnumarker=5
                if countH2 == 1 and countO2 == 1:
                    fg_key='F-CH-O'
                    ldalabel=21
                    gnucolor=14423100
                    gnumarker=7
                #aldehyde
                if countC2 == 1 and countO2 == 1:
                    fg_key='F-C-CO'
                    ldalabel=22
                    gnucolor=2003199
                    gnumarker=13
                #ketone
                if countC2 == 1 and countO2 == 1:
                    fg_key='F-C-CO'
                    ldalabel=23
                    gnucolor=3100495
                    gnumarker=5
                #ester/carboxylicacid
                if countO2 == 2:
                    fg_key='F-C-OO'
                    ldalabel=24
                    gnucolor=6737322
                    gnumarker=3
                if countN2 == 1 and countO2 ==1:
                    fg_key='F-C-NO'
                    ldalabel=26
                    gnucolor=16776960
                    gnumarker=9
                if countC2 == 1 and countN2 == 1:
                    fg_key='F-C-CN'
                    ldalabel=27
                    gnucolor=16766720
                    gnumarker=1
                if countN2 == 2:
                    fg_key='F-C-NN'
                    ldalabel=28
                    gnucolor=16711935
                    gnumarker=3
                if countH2 == 1 and countN2 == 1:
                    fg_key='F-CH-N'
                    ldalabel=29
                    gnucolor=16738740
                    gnumarker=7
                if countC2 == 1 and countF2 == 2:
                    fg_key='F-C-CF'
                    ldalabel=36
                    gnucolor=16761035
                    gnumarker=3
                if  countN2 == 1 and countF2 == 2:
                    fg_key='F-C-NF'
                    ldalabel=37
                    gnucolor=12357519
                    gnumarker=13
                ##
        #                    if countH == 3:
        #                        label = 47
                
                
            if total2 == 2:
                if countC2 == 1:
                    fg_key='F-C-C'
                    ldalabel=39
                    gnucolor=8421376
                    gnumarker=7
                if countN2 == 1:
                    fg_key='F-C-N'
                    ldalabel=41
                    gnucolor=14596231
                    gnumarker=5


    if ldalabel == 999:
        print('999')
        error_message = 'this F-type functional group is unknown, molecule_id = %s atom_index =  %s, H, %s, C, %s, N, %s, O, %s, F, %s' %(each_molecule,each_atom,countH,countC,countN,countO,countF)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker