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
            fg_key = 'OH-C'
            ldalabel=3
            gnucolor=9498256
            gnumarker=3
            hexcolor='#90EE90'
        if countH == 1 and countN == 1:
            fg_key='OH-N'
            ldalabel = 4
            gnucolor=15787660
            gnumarker=7 
            hexcolor='#F0E68C' 
        if countN == 2:
            fg_key='O-NN'
            ldalabel = 5
            gnucolor=16766720
            gnumarker=5
            hexcolor='#FFD700' 
    if total == 1:
        if countC == 1:
            fg_key='O-C'
            ldalabel=6
            gnucolor=14315734
            gnumarker=2
            hexcolor='#DA70D6'
        if countN == 1:
            fg_key='O-N'
            ldalabel= 7
            gnucolor=10824234
            gnumarker=7
            hexcolor='#A52A2A' 

    if ldalabel == 999:
        print('999')
        error_message = 'this O-type functional group is unknown, molecule_id = %s atom_index =  %s, H, %s %s, C, %s %s, N, %s %s, O, %s %s, F, %s %s' %(each_molecule,each_atom,countH,countH2,countC,countC2,countN,countN2,countO,countO2,countF,countF2)
        raise ValueError(error_message)
            
    return fg_key,ldalabel,gnucolor,gnumarker,hexcolor