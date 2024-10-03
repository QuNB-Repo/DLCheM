import random
import os
import fileinput
import pandas as pd
import copy
import re   
    
def crossover(position_N_to_add,position_db_to_add,position_fg_to_add,fg_to_add,linker_to_add)
    parent1, parent2 = random.choices(os.listdir('./bpy_GA/mol-files/'), k=2)
    parent1 = str(parent1)
    parent2 = str(parent2)

    #SECTIONS OF PARENT 1
    start = parent1.index('N')
    end = parent1.index('D',start+1)
    N_1 = parent1[start+1:end]

    start = parent1.index('B')
    end = parent1.index('F',start+1)
    DB_1 = parent1[start+1:end]

    start = parent1.index('G')
    end = parent1.index(',',start+1)
    FG_pos_1 = parent1[start+1:end]

    start = parent1.index(',')
    end = parent1.index('L',start+1)
    FG_1 = parent1[start+1:end]

    start = parent1.index('L')
    end = parent1.index('.',start+1)
    L_1 = parent1[start+1:end]

    #SECTIONS OF PARENT 2
    start = parent2.index('N')
    end = parent2.index('D',start+1)
    N_2 = parent2[start+1:end]

    start = parent2.index('B')
    end = parent2.index('F',start+1)
    DB_2 = parent2[start+1:end]

    start = parent2.index('G')
    end = parent2.index(',',start+1)
    FG_pos_2 = parent2[start+1:end]

    start = parent2.index(',')
    end = parent2.index('L',start+1)
    FG_2 = parent2[start+1:end]

    start = parent2.index('L')
    end = parent2.index('.',start+1)
    L_2 = parent2[start+1:end]

    def crossover_uniform(gene1,gene2,prop):
        c1 = copy.deepcopy(gene1)
        c2 = copy.deepcopy(gene2)

        for i in range(len(gene1)):
            if random.random() < prop:
                c1[i], c2[i] = gene2[i], gene1[i]

        return[c1, c2]

    random.seed()

    ##[N atom position, double bond position, substituent position, substituent, linker]
    gene1 = ['N'+N_1,'DB'+DB_1,'FG'+FG_pos_1,','+FG_1,'L'+L_1]
    gene2 = ['N'+N_2,'DB'+DB_2,'FG'+FG_pos_2,','+FG_2,'L'+L_2]

    offspring = crossover_uniform(gene1, gene2, 0.5)
    child1 = offspring[0]
    child2 = offspring[1]

    print('Parent 1:',parent1)
    print('Parent 2:',parent2)
    print('Child 1:',child1[0]+child1[1]+child1[2]+child1[3]+child1[4])
    print('Child 2:',child2[0]+child2[1]+child2[2]+child2[3]+child2[4])
        
    #N atom positions
    position_N_to_add = child1[0]
    position_N_to_add = str(position_N_to_add)
    position_N_to_add = re.findall('[0-9]+',position_N_to_add)
    position_N_to_add = [eval(i) for i in position_N_to_add]

    #double bond positions
    position_db_to_add = child1[1]
    position_db_to_add = str(position_db_to_add)
    position_db_to_add = re.findall('[0-9]+',position_db_to_add)
    position_db_to_add = [eval(i) for i in position_db_to_add]

    #fgroup positions
    position_fg_to_add = child1[2]
    position_fg_to_add = str(position_fg_to_add)
    position_fg_to_add = re.findall('[0-9]+',position_fg_to_add)
    position_fg_to_add = [eval(i) for i in position_fg_to_add]

    #fgroups
    fg_to_add = child1[3]
    fg_to_add = str(fg_to_add)
    fg_to_add = re.findall('[0-9]+',fg_to_add)
    fg_to_add = [eval(i) for i in fg_to_add]

    #linker
    linker_to_add = child1[4]
    linker_to_add = str(linker_to_add)
    linker_to_add = re.findall('[0-9]+',linker_to_add)
    linker_to_add = [eval(i) for i in linker_to_add]

    if linker_to_add == 100:
        double = True
    else:
        double = False
        
    print('N atoms:', position_N_to_add)
    print('double bonds:', position_db_to_add)
    print('fgroup position:', position_fg_to_add)
    print('fgroup:', fg_to_add)
    print('linker', linker_to_add)
    
    return position_N_to_add,position_db_to_add,position_fg_to_add,fg_to_add,linker_to_add