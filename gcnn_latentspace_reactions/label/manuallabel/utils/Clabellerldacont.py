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
            ldalabel=54
            gnucolor=2263842 
            gnumarker=13    
        if countH == 1 and countC == 3:
            fg_key='CH-CCC'
            ldalabel=55
            gnucolor=10025880
            gnumarker=7
        if countH == 2 and countC == 2:
            fg_key='CHH-CC'
            ldalabel=56
            gnucolor=10145074
            gnumarker=5
        if countH == 3 and countC == 1:
            fg_key='CHHH-C'
            ldalabel=57
            gnucolor=3329330
            gnumarker=1
        if countH == 4: 
            fg_key='CHHHH'
            ldalabel=58
            gnucolor=16711680
            gnumarker=9
        if countH == 2 and countC == 1 and countO == 1:
            fg_key='CHH-CO'
            ldalabel=59
            gnucolor=11529966
            gnumarker=1
        if countH == 2 and countC == 1 and countN == 1:
            fg_key='CHH-CN'
            ldalabel=60
            gnucolor=9109643
            gnumarker=3
        if countH == 3 and countO == 1:
            fg_key='CHHH-O'
            ldalabel=61
            gnucolor=4620980
            gnumarker=2
        if countH == 3 and countN == 1:
            fg_key='CHHH-N'
            ldalabel=62
            gnucolor=8388736
            gnumarker=7
        if countH == 1 and countC == 2 and countO == 1:
            fg_key='CH-CCO'
            ldalabel=63
            gnucolor=35723
            gnumarker=11
        if countH == 1 and countC == 1 and countO == 2:
            fg_key='CH-COO'
            ldalabel=64
            gnucolor=6591981
            gnumarker=2
        if countH == 1 and countC == 2 and countN == 1:
            fg_key='CH-CCN'
            ldalabel=65
            gnucolor=1644912
            gnumarker=3
        if countC == 3 and countO == 1:
            fg_key='C-CCCO'
            ldalabel=66
            gnucolor=8900331
            gnumarker=15
        if countC == 3 and countN == 1:
            fg_key='C-CCCN'
            ldalabel=67
            gnucolor=4734347
            gnumarker=11
        if countH == 2 and countO == 2:
            fg_key='CHH-OO'
            ldalabel=68
            gnucolor=4286945
            gnumarker=5
        if countF == 4: 
            fg_key='CF4'
            ldalabel=69
            gnucolor=11674146
            gnumarker=15
        if countF == 3 and countC == 1:
            fg_key='CF3-C'
            ldalabel=70
            gnucolor=13047173
            gnumarker=1
        if countC == 2 and countO == 2:
            fg_key='C-CCOO'
            ldalabel=71
            gnucolor=2142890
            gnumarker=3
    if total == 3:
        if countC == 3:
            fg_key='C-CCC'
            ldalabel=72
            gnucolor=9419919
            gnumarker=13
        if countH == 1 and countC == 2:
            fg_key='CH-CC'
            ldalabel=73
            gnucolor=32768
            gnumarker=5
        if countH == 2 and countO == 1:
            fg_key='CHH-O'
            ldalabel=74
            gnucolor=14423100
            gnumarker=7
        #aldehyde
        if countH == 1 and countC == 1 and countO == 1:
            fg_key='CH-CO'
            ldalabel=75
            gnucolor=2003199
            gnumarker=13
        #ketone
        if countC == 2 and countO == 1:
            fg_key='C-CCO'
            ldalabel=76
            gnucolor=3100495
            gnumarker=5
        #ester/carboxylicacid
        if countC == 1 and countO == 2:
            fg_key='C-COO'
            ldalabel=77
            gnucolor=6737322
            gnumarker=3
        if countH == 1 and countO == 2:
            fg_key='CH-OO'
            ldalabel=78
            gnucolor=11584734
            gnumarker=5
        if countH == 1 and countN == 1 and countO ==1:
            fg_key='CH-NO'
            ldalabel=79
            gnucolor=16776960
            gnumarker=9
        if countC == 1 and countN == 1 and countO ==1:
            fg_key='C-CNO'
            ldalabel=80
            gnucolor=16766720
            gnumarker=1
        if countC == 1 and countN == 2:
            fg_key='C-CNN'
            ldalabel=81
            gnucolor=16711935
            gnumarker=3
        if countH == 1 and countC == 1 and countN == 1:
            fg_key='CH-CN'
            ldalabel=82
            gnucolor=16738740
            gnumarker=7
        if countN == 2 and countO == 1:
            fg_key='C-NNO'
            ldalabel=83
            gnucolor=16716947
            gnumarker=11
        if countN == 1 and countO == 2:
            fg_key='C-NOO'
            ldalabel=84
            gnucolor=14329120
            gnumarker=15
        if countN == 3:
            fg_key='C-NNN'
            ldalabel=85
            gnucolor=9662683
            gnumarker=2
        if countO == 3: 
            fg_key='C-OOO'
            ldalabel=86
            gnucolor=16753920
            gnumarker=13
        if countC == 2 and countN == 1:
            fg_key='C-CCN'
            ldalabel=87
            gnucolor=8087790
            gnumarker=1
        if countH == 1 and countN == 2:
            fg_key='CH-NN'
            ldalabel=88
            gnucolor=14381203
            gnumarker=15
        if countC == 2 and countF == 1:
            fg_key='C-CCF'
            ldalabel=89
            gnucolor=16761035
            gnumarker=3
        if countC == 1 and countN == 1 and countF == 1:
            fg_key='C-CNF'
            ldalabel=90
            gnucolor=12357519
            gnumarker=13
        if countN == 2 and countF == 1:
            fg_key='C-NNF'
            ldalabel=91
            gnucolor=13458524
            gnumarker=15

        ##
#                    if countH == 3:
#                        label = 47
        #NEWWW!!! COLOR NOT ASSIGNED
#        if countH == 2 and countC == 1:
#            fg_key = 'CHH-C'
#            ldalabel = 43
#            gnucolor = 000000
#            gnumarker = 1

        
        
    if total == 2:
        if countC == 2:
            fg_key='C-CC'
            ldalabel=92
            gnucolor=8421376
            gnumarker=7
        if countH == 1 and countC == 1:
            fg_key='CH-C'
            ldalabel=93
            gnucolor=16729344
            gnumarker=7
        if countH == 1 and countN == 1:
            fg_key='CH-N'
            ldalabel=94
            gnucolor=14596231
            gnumarker=5
        if countC == 1 and countN == 1:
            fg_key='C-CN'
            ldalabel=95
            gnucolor=13468991
            gnumarker=3

        #THIS IS NEWWWW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        #AND MAY BE A MISTAKE FROM RDKIT 
#        if countN == 2:
#            fg_key = 'C-NN'
#            ldalabel=42
#            gnucolor=13789470
#            gnumarker=3
        ##
#        if countO == 2:
#            label = 42
        ##
#        if countC ==1 and countO == 1:
#            label = 43
    
    #these are extra
#    if total == 1:
        ##
#        if countO == 1:
#            label = 44
        ##
#        if countC == 1:
#            label=45
        #
#        if countN == 1:
#            label = 46

    if ldalabel == 999:
        print('999')
        error_message = 'this C-type functional group is unknown, molecule_id = %s, atom_index = %s, H, %s, C, %s, N, %s, O, %s, F, %s' %(each_molecule,each_atom,countH,countC,countN,countO,countF)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker
