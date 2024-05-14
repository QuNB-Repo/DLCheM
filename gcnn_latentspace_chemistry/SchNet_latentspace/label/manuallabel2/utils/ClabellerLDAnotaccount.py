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
                if '9' in each_line[0:1]:
                    countF = countF+1
                    
    total = countC + countH + countO + countN + countF

    ldalabel=999

    if total == 4:
        if countC == 4:
            fg_key='C-CCCC'
            ldalabel=0
            gnucolor=2263842 
            gnumarker=13    
        if countH == 1 and countC == 3:
            fg_key='CH-CCC'
            ldalabel=1
            gnucolor=10025880
            gnumarker=7
        if countH == 2 and countC == 2:
            fg_key='CHH-CC'
            ldalabel=2
            gnucolor=10145074
            gnumarker=5
        if countH == 3 and countC == 1:
            fg_key='CHHH-C'
            ldalabel=3
            gnucolor=3329330
            gnumarker=1
        if countH == 4: 
            fg_key='CHHHH'
            ldalabel=4
            gnucolor=16711680
            gnumarker=9
        if countH == 2 and countC == 1 and countO == 1:
            fg_key='CHH-CO'
            ldalabel=5
            gnucolor=11529966
            gnumarker=1
        if countH == 2 and countC == 1 and countN == 1:
            fg_key='CHH-CN'
            ldalabel=6
            gnucolor=9109643
            gnumarker=3
        if countH == 3 and countO == 1:
            fg_key='CHHH-O'
            ldalabel=7
            gnucolor=4620980
            gnumarker=2
        if countH == 3 and countN == 1:
            fg_key='CHHH-N'
            ldalabel=8
            gnucolor=8388736
            gnumarker=7
        if countH == 1 and countC == 2 and countO == 1:
            fg_key='CH-CCO'
            ldalabel=9
            gnucolor=35723
            gnumarker=11
        if countH == 1 and countC == 1 and countO == 2:
            fg_key='CH-COO'
            ldalabel=10
            gnucolor=6591981
            gnumarker=2
        if countH == 1 and countC == 2 and countN == 1:
            fg_key='CH-CCN'
            ldalabel=11
            gnucolor=1644912
            gnumarker=3
        if countC == 3 and countO == 1:
            fg_key='C-CCCO'
            ldalabel=12
            gnucolor=8900331
            gnumarker=15
        if countC == 3 and countN == 1:
            fg_key='C-CCCN'
            ldalabel=13
            gnucolor=4734347
            gnumarker=11
        if countH == 2 and countO == 2:
            fg_key='CHH-OO'
            ldalabel=14
            gnucolor=4286945
            gnumarker=5
        if countF == 4: 
            fg_key='CF4'
            ldalabel=15
            gnucolor=11674146
            gnumarker=15
        if countF == 3 and countC == 1:
            fg_key='CF3-C'
            ldalabel=16
            gnucolor=13047173
            gnumarker=1
        if countC == 2 and countO == 2:
            fg_key='C-CCOO'
            ldalabel=17
            gnucolor=2142890
            gnumarker=3
    if total == 3:
        if countC == 3:
            fg_key='C-CCC'
            ldalabel=18
            gnucolor=9419919
            gnumarker=13
        if countH == 1 and countC == 2:
            fg_key='CH-CC'
            ldalabel=19
            gnucolor=32768
            gnumarker=5
        if countH == 2 and countO == 1:
            fg_key='CHH-O'
            ldalabel=20
            gnucolor=14423100
            gnumarker=7
        #aldehyde
        if countH == 1 and countC == 1 and countO == 1:
            fg_key='CH-CO'
            ldalabel=21
            gnucolor=2003199
            gnumarker=13
        #ketone
        if countC == 2 and countO == 1:
            fg_key='C-CCO'
            ldalabel=22
            gnucolor=3100495
            gnumarker=5
        #ester/carboxylicacid
        if countC == 1 and countO == 2:
            fg_key='C-COO'
            ldalabel=23
            gnucolor=6737322
            gnumarker=3
        if countH == 1 and countO == 2:
            fg_key='CH-OO'
            ldalabel=24
            gnucolor=11584734
            gnumarker=5
        if countH == 1 and countN == 1 and countO ==1:
            fg_key='CH-NO'
            ldalabel=25
            gnucolor=16776960
            gnumarker=9
        if countC == 1 and countN == 1 and countO ==1:
            fg_key='C-CNO'
            ldalabel=26
            gnucolor=16766720
            gnumarker=1
        if countC == 1 and countN == 2:
            fg_key='C-CNN'
            ldalabel=27
            gnucolor=16711935
            gnumarker=3
        if countH == 1 and countC == 1 and countN == 1:
            fg_key='CH-CN'
            ldalabel=28
            gnucolor=16738740
            gnumarker=7
        if countN == 2 and countO == 1:
            fg_key='C-NNO'
            ldalabel=29
            gnucolor=16716947
            gnumarker=11
        if countN == 1 and countO == 2:
            fg_key='C-NOO'
            ldalabel=30
            gnucolor=14329120
            gnumarker=15
        if countN == 3:
            fg_key='C-NNN'
            ldalabel=31
            gnucolor=9662683
            gnumarker=2
        if countO == 3: 
            fg_key='C-OOO'
            ldalabel=32
            gnucolor=16753920
            gnumarker=13
        if countC == 2 and countN == 1:
            fg_key='C-CCN'
            ldalabel=33
            gnucolor=8087790
            gnumarker=1
        if countH == 1 and countN == 2:
            fg_key='CH-NN'
            ldalabel=34
            gnucolor=14381203
            gnumarker=15
        if countC == 2 and countF == 1:
            fg_key='C-CCF'
            ldalabel=35
            gnucolor=16761035
            gnumarker=3
        if countC == 1 and countN == 1 and countF == 1:
            fg_key='C-CNF'
            ldalabel=36
            gnucolor=12357519
            gnumarker=13
        if countN == 2 and countF == 1:
            fg_key='C-NNF'
            ldalabel=37
            gnucolor=13458524
            gnumarker=15

    if total == 2:
        if countC == 2:
            fg_key='C-CC'
            ldalabel=38
            gnucolor=8421376
            gnumarker=7
        if countH == 1 and countC == 1:
            fg_key='CH-C'
            ldalabel=39
            gnucolor=16729344
            gnumarker=7
        if countH == 1 and countN == 1:
            fg_key='CH-N'
            ldalabel=40
            gnucolor=14596231
            gnumarker=5
        if countC == 1 and countN == 1:
            fg_key='C-CN'
            ldalabel=41
            gnucolor=13468991
            gnumarker=3

    #THESE ARE THE MOST RECENT ADDITIONS!!! NOT UPDATED IN LEGEND! 
    if total == 4:
        if countH == 1 and countC == 1 and countN == 2:
            fg_key = 'CH-CNN'
            ldalabel = 43
            gnucolor = 16772045
            gnumarker = 11
        if countH == 1 and countC == 1 and countN == 1 and countO == 1:
            fg_key = 'CH-CNO'
            ldalabel = 44
            gnucolor = 16119260
            gnumarker = 2
        if countC ==  1 and countO == 3:
            fg_key = 'C-COOO'
            ldalabel = 45
            gnucolor = 13789470
            gnumarker = 7
        if countC == 2 and countN == 1 and countO == 1:
            fg_key = 'C-CCNO'
            ldalabel = 46
            gnucolor = 8620085
            gnumarker = 7
    if total == 3:
        ##NEWWW: FOUND AT 21581 in QM9
        if countC == 1 and countO == 1 and countF == 1:
            fg_key = 'C-COF'
            ldalabel = 47
            gnucolor = 7048739
            gnumarker = 9
        ##NEWWW: FOUND AT 21581 in QM9
        if countN == 1 and countO == 1 and countF == 1:
            fg_key = 'C-NOF'
            ldalabel = 48
            gnucolor = 5597999
            gnumarker = 1
        if countH == 2 and countC == 1:
            fg_key = 'CH2-C'
            ldalabel = 49
            gnucolor = 16774638
            gnumarker = 5
    if total == 2:
        if countC == 1 and countO == 1:
            fg_key = 'C-CO'
            ldalabel = 50
            gnucolor = 1446290
            gnumarker = 15
        if countN == 1 and countO == 1:
            fg_key = 'C-NO'
            ldalabel = 51
            gnucolor = 16032864
            gnumarker = 3
        if countO == 2:
            fg_key = 'C-OO'
            ldalabel = 52
            gnucolor = 15794175
            gnumarker = 2
    if total == 1:
        if countN == 1:
            fg_key = 'C-N'
            ldalabel = 53
            gnucolor = 11591910
            gnumarker = 9
        if countO == 1:
            fg_key == 'C-O'
            ldalabel = 54
            gnucolor = 1136921
            gnumarker = 3



    if ldalabel == 999:
        print('999')
        error_message = 'this C-type functional group is unknown, molecule_id = %s, atom_index = %s, H, %s, C, %s, N, %s, O, %s, F, %s' %(each_molecule,each_atom,countH,countC,countN,countO,countF)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker
