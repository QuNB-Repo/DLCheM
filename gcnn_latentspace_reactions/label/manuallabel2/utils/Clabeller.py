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
            hexcolor='#228B22'
            gnumarker=13    
        if countH == 1 and countC == 3:
            fg_key='CH-CCC'
            ldalabel=1
            gnucolor=10025880
            gnumarker=7
            hexcolor='#98FB98'
        if countH == 2 and countC == 2:
            fg_key='CHH-CC'
            ldalabel=2
            gnucolor=10145074
            gnumarker=5
            hexcolor='#9ACD32'
        if countH == 3 and countC == 1:
            fg_key='CHHH-C'
            ldalabel=3
            gnucolor=3329330
            gnumarker=1
            hexcolor='#32CD32'
        if countH == 4: 
            fg_key='CHHHH'
            ldalabel=4
            gnucolor=16711680
            gnumarker=9
            hexcolor='#FF0000'
        if countH == 2 and countC == 1 and countO == 1:
            fg_key='CHH-CO'
            ldalabel=5
            gnucolor=11529966
            gnumarker=1
            hexcolor='#AFEEEE'
        if countH == 2 and countC == 1 and countN == 1:
            fg_key='CHH-CN'
            ldalabel=6
            gnucolor=9109643
            gnumarker=3
            hexcolor='#8B008B'
        if countH == 3 and countO == 1:
            fg_key='CHHH-O'
            ldalabel=7
            gnucolor=4620980
            gnumarker=2
            hexcolor='#4682B4'
        if countH == 3 and countN == 1:
            fg_key='CHHH-N'
            ldalabel=8
            gnucolor=8388736
            gnumarker=7
            hexcolor='#800080'
        if countH == 1 and countC == 2 and countO == 1:
            fg_key='CH-CCO'
            ldalabel=9
            gnucolor=35723
            gnumarker=11
            hexcolor='#008B8B'
        if countH == 1 and countC == 1 and countO == 2:
            fg_key='CH-COO'
            ldalabel=10
            gnucolor=6591981
            gnumarker=2
            hexcolor='#6495ED'
        if countH == 1 and countC == 2 and countN == 1:
            fg_key='CH-CCN'
            ldalabel=11
            gnucolor=1644912
            gnumarker=3
            hexcolor='#191970'
        if countC == 3 and countO == 1:
            fg_key='C-CCCO'
            ldalabel=12
            gnucolor=8900331
            gnumarker=15
            hexcolor='#87CEEB'
        if countC == 3 and countN == 1:
            fg_key='C-CCCN'
            ldalabel=13
            gnucolor=4734347
            gnumarker=11
            hexcolor='#483D8B'
        if countH == 2 and countO == 2:
            fg_key='CHH-OO'
            ldalabel=14
            gnucolor=4286945
            gnumarker=5
            hexcolor='#4169E1'
        if countF == 4: 
            fg_key='CF4'
            ldalabel=15
            gnucolor=11674146
            gnumarker=15
            hexcolor='#B22222'
        if countF == 3 and countC == 1:
            fg_key='CF3-C'
            ldalabel=16
            gnucolor=13047173
            gnumarker=1
            hexcolor='#C71585'
        if countC == 2 and countO == 2:
            fg_key='C-CCOO'
            ldalabel=17
            gnucolor=2142890
            gnumarker=3
            hexcolor='#20B2AA'
    if total == 3:
        if countC == 3:
            fg_key='C-CCC'
            ldalabel=18
            gnucolor=9419919
            gnumarker=13
            hexcolor='#8FBC8F'
        if countH == 1 and countC == 2:
            fg_key='CH-CC'
            ldalabel=19
            gnucolor=32768
            gnumarker=5
            hexcolor='#008000'
        if countH == 2 and countO == 1:
            fg_key='CHH-O'
            ldalabel=20
            gnucolor=14423100
            gnumarker=7
            hexcolor='#DC143C'
        #aldehyde
        if countH == 1 and countC == 1 and countO == 1:
            fg_key='CH-CO'
            ldalabel=21
            gnucolor=2003199
            gnumarker=13
            hexcolor='#1E90FF'
        #ketone
        if countC == 2 and countO == 1:
            fg_key='C-CCO'
            ldalabel=22
            gnucolor=3100495
            gnumarker=5
            hexcolor='#2F4F4F'
        #ester/carboxylicacid
        if countC == 1 and countO == 2:
            fg_key='C-COO'
            ldalabel=23
            gnucolor=6737322
            hexcolor='#66CDAA'
            gnumarker=3
        if countH == 1 and countO == 2:
            fg_key='CH-OO'
            ldalabel=24
            gnucolor=11584734
            hexcolor='#B0C4DE'
            gnumarker=5
        if countH == 1 and countN == 1 and countO ==1:
            fg_key='CH-NO'
            ldalabel=25
            gnucolor=16776960
            gnumarker=9
            hexcolor='#FFFF00'
        if countC == 1 and countN == 1 and countO ==1:
            fg_key='C-CNO'
            ldalabel=26
            gnucolor=16766720
            gnumarker=1
            hexcolor='#FFD700'
        if countC == 1 and countN == 2:
            fg_key='C-CNN'
            ldalabel=27
            gnucolor=16711935
            gnumarker=3
            hexcolor='#FF00FF'
        if countH == 1 and countC == 1 and countN == 1:
            fg_key='CH-CN'
            ldalabel=28
            gnucolor=16738740
            gnumarker=7
            hexcolor='#FF69B4'
        if countN == 2 and countO == 1:
            fg_key='C-NNO'
            ldalabel=29
            gnucolor=16716947
            gnumarker=11
            hexcolor='#FF1493'
        if countN == 1 and countO == 2:
            fg_key='C-NOO'
            ldalabel=30
            gnucolor=14329120
            gnumarker=15
            hexcolor='#DAA520'
        if countN == 3:
            fg_key='C-NNN'
            ldalabel=31
            gnucolor=9662683
            gnumarker=2
            hexcolor='#9370DB'
        if countO == 3: 
            fg_key='C-OOO'
            ldalabel=32
            gnucolor=16753920
            gnumarker=13
            hexcolor='#FFA500'
        if countC == 2 and countN == 1:
            fg_key='C-CCN'
            ldalabel=33
            gnucolor=8087790
            gnumarker=1
            hexcolor='#7B68EE'
        if countH == 1 and countN == 2:
            fg_key='CH-NN'
            ldalabel=34
            gnucolor=14381203
            gnumarker=15
            hexcolor='#DB7093'
        if countC == 2 and countF == 1:
            fg_key='C-CCF'
            ldalabel=35
            gnucolor=16761035
            gnumarker=3
            hexcolor='#FFC0CB'
        if countC == 1 and countN == 1 and countF == 1:
            fg_key='C-CNF'
            ldalabel=36
            gnucolor=12357519
            gnumarker=13
            hexcolor='#BC8F8F'
        if countN == 2 and countF == 1:
            fg_key='C-NNF'
            ldalabel=37
            gnucolor=13458524
            gnumarker=15
            hexcolor='#CD5C5C'

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
            ldalabel=38
            gnucolor=8421376
            gnumarker=7
            hexcolor='#808000'
        if countH == 1 and countC == 1:
            fg_key='CH-C'
            ldalabel=39
            gnucolor=16729344
            gnumarker=7
            hexcolor='#FF4500'
        if countH == 1 and countN == 1:
            fg_key='CH-N'
            ldalabel=40
            gnucolor=14596231
            gnumarker=5
            hexcolor='#DEB887'
        if countC == 1 and countN == 1:
            fg_key='C-CN'
            ldalabel=41
            gnucolor=13468991
            gnumarker=3
            hexcolor='#CD853F'

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
            
    return fg_key,ldalabel,gnucolor,gnumarker,hexcolor
