    

from label.manuallabel2.utils import utils


#label function for H elements
def label(connections_to_atom,xyz_file_read,mol_file_read,number_atoms,each_molecule,each_atom):

    #count the atoms around the hydrogen, usually there is just one! 
    countC = 0
    countH = 0
    countO = 0
    countN = 0
    countF = 0
    countH2 = 0
    countC2 = 0
    countN2 = 0
    countO2 = 0
    countF2 = 0
    countH3 = 0
    countC3 = 0 
    countN3 = 0
    countO3 = 0 
    countF3 = 0
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

    #figure out if it is C, N, or O
    C_nbh = False
    N_nbh = False
    O_nbh = False
    
    if countC == 1:
        C_nbh = True
    if countO == 1:
        O_nbh = True
    if countN == 1:
        N_nbh = True

    #usually there is just 1 connection to hydrogen, so neighbor is just min of a single array
    neighbor_of_hydrogen = int(min(connections_to_atom))

    #get the neighbors of carbon/nitrogen/oxygen that next to the hyrogen
    neighbors_of_CNO = utils.find_connections(number_atoms,neighbor_of_hydrogen,mol_file_read)
    total2,countC2,countH2,countO2,countN2, countF2 = utils.count_nn(neighbors_of_CNO,xyz_file_read)
    
    ldalabel=999
    if C_nbh == True:
        if total2 == 4:    
            if countH2 == 4:
                #this is methane, no cluster
                fg_key='H-C-HHH'
                ldalabel=0
                gnucolor=16711680 
                gnumarker=13     
            elif countH2==3 and countC2==1:
                #these are all the methyl groups
                fg_key = 'H-CHH-C'
                ldalabel=1
                gnucolor = 10025880  
                gnumarker = 7                   
            #methyl-O
            elif countH2 == 3 and countO2 == 1:
                fg_key = 'H-CHH-O'
                ldalabel=2
                gnucolor = 10145074  
                gnumarker = 5         
            #methyl-N
            elif countH2 == 3 and countN2 == 1:
                fg_key = 'H-CHH-N'
                ldalabel=3
                gnucolor = 3329330    
                gnumarker = 1                     
            #this is methylene
            elif countH2 == 2 and countC2 == 2:
                fg_key = 'H-CH-CC'
                ldalabel=4
                gnucolor = 11529966    
                gnumarker = 9                  
            #methylene-N
            elif countH2==2 and countC2==1 and countN2==1:
                fg_key = 'H-CH-CN'
                ldalabel=5
                gnucolor = 3978097    
                gnumarker = 1   
            elif countH2 == 2 and countC2 == 1 and countO2==1:
                fg_key = 'H-CH-CO'
                ldalabel = 6
                gnucolor = 9109643    
                gnumarker = 3                         
            elif countH2 == 2 and countO2 == 2:
                fg_key = 'H-CH-OO'
                ldalabel = 7
                gnucolor = 4620980    
                gnumarker = 2         
            #methines
            elif countH2==1 and countC2==3:
                fg_key='H-C-CCC'
                ldalabel=8
                gnucolor=3100495      
                gnumarker=7         
            elif countH2==1 and countC2==2 and countN2==1:
                fg_key='H-C-CCN'
                ldalabel= 9
                gnucolor= 35723
                gnumarker= 11 
            elif countH2==1 and countC2==2 and countO2==1:
                fg_key = 'H-C-CCO'
                ldalabel=10
                gnucolor= 4734347 
                gnumarker=2  
            #outer acetal
            elif countH2==1 and countC2==1 and countO2==2:
                fg_key = 'H-CH-COO'
                ldalabel=11 
                gnucolor=1644912       
                gnumarker=3           
        elif total2 == 3:
            #trigonal planar carbons
            #arromatic
            if countH2==1 and countC2==2:
                fg_key = 'H-C-CC'
                ldalabel=12
                gnucolor=9699539      
                gnumarker=15          
            elif countH2==1 and countC2==1 and countO2==1:
                fg_key = 'H-C-CO'
                ldalabel=13
                gnucolor=8421376    
                gnumarker=11         
            elif countH2==1 and countC2==1 and countN2==1:
                fg_key = 'H-C-CN'
                ldalabel=14
                gnucolor=4169E1    
                gnumarker=5   
            elif countH2==1 and countO2==2:
                fg_key = 'H-C-OO'
                ldalabel=15
                gnucolor=14315734   
                gnumarker=15   
            elif countH2==1 and countN2==2:
                fg_key = 'H-C-NN'
                ldalabel=16
                gnucolor=9055202 
                gnumarker=1
            elif countH2==1 and countN2==1 and countO2==1:
                fg_key = 'H-C-NO'
                ldalabel=17
                gnucolor=6266528
                gnumarker=3
            elif countH2==2 and countO2==1:
                fg_key = 'H-CH-O'
                ldalabel=18
                gnucolor=11674146
                gnumarker=13
                #formaldehyde
        elif total2 == 2:
            if countH2==1 and countC2==1:
                fg_key = 'H-C-C'
                ldalabel=19
                gnucolor=11674146
                gnumarker=13
            elif countH2==1 and countN2==1:
                fg_key = 'H-C-N'
                ldalabel=20
                gnucolor=14204888
                gnumarker=7               
    if N_nbh == True:

        if total2 == 3:
            #trigonal planar nitrogens
            #first one is ammonia, no cluster
            #ammonia            
            if countH2==3:
                fg_key = 'H-NHH'
                ldalabel=21
                gnucolor=14423100
                gnumarker=2     
            #primary amine
            elif countH2==2 and countC2==1:

                #locate the carbon neighbor
                check_if_hydrogen = True
                while check_if_hydrogen == True:
                    check_neighbor = min(neighbors_of_CNO)
                    element_present = utils.check_element(int(check_neighbor),xyz_file_read)

                    if element_present == 'H':
                        check_if_hydrogen = True
                        index = neighbors_of_CNO.index(int(check_neighbor))
                        neighbors_of_CNO.pop(index)

                    else:
                        carbon_neighbor = check_neighbor
                        check_if_hydrogen = False

                #get the neighbors of this carbon
                carbon_neighbor_neighbors = utils.find_connections(number_atoms,carbon_neighbor,mol_file_read)
                total3,countC3,countH3,countO3,countN3,countF3 = utils.count_nn(carbon_neighbor_neighbors,xyz_file_read)

                if total3 == 4:
                    if countC3 == 3:
                        fg_key = 'H-NH-C-CCC'
                        ldalabel=22
                        gnucolor=15657130
                        gnumarker=5
                    if countC3 == 2 and countH3 == 1:
                        fg_key = 'H-NH-CH-CC'
                        ldalabel=23
                        gnucolor=3468991
                        gnumarker=3
                    if countC3 == 1 and countH3 == 2:
                        fg_key = 'H-NH-CHH-C'
                        ldalabel=24
                        gnucolor=16768685
                        gnumarker=5


                if total3 == 3:
                    if countC3 == 2: 
                        fg_key = 'H-NH-C-CC'
                        ldalabel=25
                        gnucolor=16776960
                        gnumarker=9
                    if countC3 == 1 and countN3 == 2:
                        fg_key = 'H-NH-C-CN'
                        ldalabel=26
                        gnucolor=16766720
                        gnumarker=1
                    if countC3 == 1 and countO3 == 1:
                        fg_key = 'H-NH-C-CN'
                        ldalabel=27
                        gnucolor=14329120
                        gnumarker=3
                    if countO3 == 2:
                        fg_key = 'H-NH-C-OO'
                        ldalabel=28
                        gnucolor=12433259
                        gnumarker=7
                    if countO3 == 1 and countH3 == 1:
                        fg_key = 'H-NH-CH-O'
                        ldalabel=29
                        gnucolor=16737095
                        gnumarker=11
                    if countO3 == 1 and countN3 == 2:
                        fg_key = 'H-NH-C-NO'
                        ldalabel=30
                        gnucolor=12092939
                        gnumarker=15
                    if countN3 == 3:
                        fg_key = 'H-NH-C-NN'
                        ldalabel=31
                        gnucolor=16032864
                        gnumarker=2

            #secondary amine
            elif countH2==1 and countC2==2:
                fg_key = 'H-N-CC'
                ldalabel=32
                gnucolor=16753920
                gnumarker=13
            elif countH2==1 and countN2==1 and countC2==1:
                fg_key = 'H-N-CN'
                ldalabel=33
                gnucolor=16747520
                gnumarker=1
            elif countH2==1 and countN2==2:
                fg_key = 'H-N-NN'
                ldalabel=34
                gnucolor=13789470
                gnumarker=15
        elif total2 == 2:
            if countH2==1 and countC2==1:
                fg_key = 'H-N-C'
                ldalabel=35
                gnucolor=16770229
                gnumarker=3
        elif total2 == 4:
            if countH2 == 3 and countC2 == 1:
                fg_key = 'H-NHHH-C'
                ldalabel=36
                gnucolor=8388736
                gnumarker=13
            if countH2 == 2 and countC2 == 2:
                fg_key = 'H-NHH-CC'
                ldalabel=37
                gnucolor=10824234
                gnumarker=15

                
    if O_nbh == True:
        if total2==2:
            #alcohols
            if countH2==1 and countC2==1:

                check_if_hydrogen = True
                while check_if_hydrogen == True:
                    check_neighbor = min(neighbors_of_CNO)
                    element_present = utils.check_element(int(check_neighbor),xyz_file_read)

                    if element_present == 'H':
                        check_if_hydrogen = True
                        index = neighbors_of_CNO.index(check_neighbor)
                        neighbors_of_CNO.pop(index)

                    else:
                        carbon_neighbor = check_neighbor
                        check_if_hydrogen = False

                #get the neighbors of this carbon
                carbon_neighbor_neighbors = utils.find_connections(number_atoms,carbon_neighbor,mol_file_read)
                total3,countC3,countH3,countO3,countN3,countF3 = utils.count_nn(carbon_neighbor_neighbors,xyz_file_read)

                if total3 == 4:
                    if countH3 == 3:
                        fg_key = 'H-O-CHHH'
                        ldalabel=38
                        gnucolor=8421504
                        gnumarker=7
                    if countC3 == 3:
                        fg_key = 'H-O-C-CCC'
                        ldalabel=39
                        gnucolor=16767673
                        gnumarker=13
                    if countC3 == 2 and countH3 == 1:
                        fg_key = 'H-O-CH-CC'
                        ldalabel=40
                        gnucolor=16752762
                        gnumarker=5
                    if countC3 == 2 and countN3 == 1:
                        fg_key = 'H-O-C-CCN'
                        ldalabel=41
                        gnucolor=16744272
                        gnumarker=3
                    if countC3 == 2 and countO3 == 2:
                        fg_key = 'H-O-C-CCO'
                        ldalabel=42
                        gnucolor=16416882
                        gnumarker=2
                    if countC3 == 1 and countH3 == 2:
                        fg_key = 'H-O-CHH-C'
                        ldalabel=43
                        gnucolor=16758465
                        gnumarker=13
                        

                if total3 == 3:
                    if countC3 == 2: 
                        fg_key = 'H-O-C-CC'
                        ldalabel=44
                        gnucolor=16738740
                        gnumarker=7
                    if countC3 == 1 and countN3 == 1:
                        fg_key = 'H-O-C-CN'
                        ldalabel=45
                        gnucolor=16716947
                        gnumarker=2
                    if countC3 == 1 and countO3 == 2:
                        fg_key = 'H-O-C-CO'
                        ldalabel=46
                        gnucolor=14381203
                        gnumarker=1
                    if countO3 == 2 and countN3 == 1:
                        fg_key = 'H-O-C-NO'
                        ldalabel=47
                        gnucolor=16711935
                        gnumarker=11
                    if countN3 == 2:
                        fg_key = 'H-O-C-NN'
                        ldalabel=48
                        gnucolor=15631086
                        gnumarker=13


            if countH2 == 2:
                fg_key = 'H-O-H'
                ldalabel=49
                gnucolor=8388608
                gnumarker=9
            if countH2==1 and countN2==1:
                fg_key = 'H-O-N'
                ldalabel=50
                gnucolor=16770229
                gnumarker=3
    ##
    ##
    ##
    ##
    ##
    #THESE ARE THE RECENT UPDATED 
    # ---> NOT INCLUDED IN CURRENT LEGEND/KEY FIGURE
    if N_nbh == True:
        if total2 == 3:            
            if countH2==2 and countC2==1:
                if total3 == 3:                    
                    if countN3 == 2 and countH3 == 1:
                        fg_key = 'H-NH-CH-N'
                        ldalabel=51
                        gnucolor=6737322
                        gnumarker=1
            if countH2 == 1 and countC2 == 1 and countO2 == 1:
                fg_key = 'H-N-CO'
                ldalabel = 52
                gnucolor = 9662683
                gnumarker = 2
        if total2 == 4:
            if countH2 == 1 and countC2 == 3:
                fg_key = 'H-NH-CCC'
                ldalabel=53
                gnucolor=5597999
                gnumarker=5
    if C_nbh == True:
        if total2 == 4: 
            if countH2 == 1 and countC2 == 1 and countN2 == 2:
                fg_key = 'H-C-CNN'
                ldalabel = 54
                gnucolor = 10506797
                gnumarker = 11
            if countH2 == 1 and countC2 == 1 and countN2 == 1 and countO2 == 1:
                fg_key = 'H-C-CNO'
                ldalabel = 55
                gnucolor = 16121850
                gnumarker = 9
        if total2 == 3:
            if countH2 == 2 and countC == 1:
                fg_key = 'H-CH-C'
                ldalabel=56
                gnucolor = 7048739
                gnumarker = 7
    if O_nbh == True:
        if total2 == 2:
            if countH2 == 1 and countC2 == 1:
                if total3 == 4: 
                    if countH3 == 1 and countC3 == 1 and countO3 == 2:
                        fg_key = 'H-O-CH-CO'
                        ldalabel = 57
                        gnucolor = 15794160
                        gnumarker = 2
                    if countH3 == 1 and countC3 == 1 and countN3 == 1 and countO3 == 1:
                        fg_key = 'H-O-CH-CN'
                        ldalabel = 58
                        gnucolor = 16444375
                        gnumarker = 11
                if total3 == 3:
                    if countH3 == 1 and countO3 == 2:
                        fg_key = 'H-O-CH-O'
                        ldalabel = 59
                        gnucolor = 12632256
                        gnumarker = 7
    #H2!! 
    if countH == 1:
        fg_key = 'H-H'
        ldalabel = 58
        gnucolor = 14745599
        gnumarker = 1

        
                
    if ldalabel == 999:
        print('999')
        error_message = 'this H-type functional group is unknown, molecule_id = %s, atom_index = %s, H, %s, %s, %s, C, %s, %s, %s, N, %s, %s, %s, O, %s, %s, %s, F, %s, %s, %s' %(each_molecule,each_atom,countH,countH2,countH3,countC,countC2,countC3,countN,countN2,countN3,countO,countO2,countO3,countF,countF2,countF3)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker