from calendar import c
import random 
import fileinput
from re import A
import numpy as np
import os

from .utils import utils_GA

import pandas as pd
import copy
import re

class autocombinatorial_GA():
    
    def __init__(self,initial_mol_filepath='',mol_string='',max_number_n=6,max_number_fg=6,max_number_db=6,double=False,molecule_name=''):

        #read mol_string from initial file
        for line in fileinput.FileInput(initial_mol_filepath,inplace=0):
            mol_string = mol_string + line

        #if double is True, dilute how often the the double linker is taken, since there is only one
        
#        if double == True:
#            linker_to_add = 100
#            double = linker_to_add

        self.mol_string = mol_string
        self.max_number_n = max_number_n
        self.max_number_fg = max_number_fg
        self.max_number_db = max_number_db
        self.double = double
        self.molecule_name = molecule_name
                        
    def show(self):
        print('max nitrogens:', self.max_number_n)
        print('double linker', self.double)
        print('molecule file:', self.mol_string)
        print('N positions', self.n_positions)
        print('fg positions', self.fgs)
        print('db position', self.db_positions)


        self.fgs = np.array(self.fgs)
        #first you have to order each one! 
        self.n_positions = np.sort(self.n_positions)
        self.db_positions = np.sort(self.db_positions)

        if len(self.fgs) > 0:
            self.fgs = self.fgs[self.fgs[:,0].argsort()]

        molecule_id = self.molecule_name + 'N'
        for i in range(len(self.n_positions)):
            molecule_id = molecule_id + '-' + str(self.n_positions[i])
        
        molecule_id = molecule_id + '-DB'
        for i in range(len(self.db_positions)):
            molecule_id = molecule_id + '-' + str(self.db_positions[i])
        
        molecule_id = molecule_id + '-FG'
        for i in range(len(self.fgs)):
            molecule_id = molecule_id + '--' + str(self.fgs[i][0]) + ',' + str(self.fgs[i][1])  

        molecule_id = molecule_id + '-L'+str(self.linker)

        self.molecule_id = molecule_id
        save_mol = open('molecules_GA/mol-files/mol-files-GA/'+molecule_id+'.mol',mode='w')
        save_mol.write(self.mol_string)

    def crossover(self,double=False):
        parent1, parent2 = random.choices(os.listdir('./molecules_GA/mol-files/'), k=2)
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
        gene1 = ['N'+N_1,'DB'+DB_1,'FG'+FG_pos_1+','+FG_1,'L'+L_1]
        gene2 = ['N'+N_2,'DB'+DB_2,'FG'+FG_pos_2+','+FG_2,'L'+L_2]

        offspring = crossover_uniform(gene1, gene2, 0.5)
        child1 = offspring[0]
        child2 = offspring[1]

        print('Parent 1:',parent1)
        print('Parent 2:',parent2)
        print('Child 1:',child1[0]+child1[1]+child1[2]+child1[3])
        print('Child 2:',child2[0]+child2[1]+child2[2]+child2[3])
        
        start = child1[2].index('G')
        end = child1[2].index(',',start+1)
        FG_pos_1 = child1[2][start+1:end]

        start = child1[2].index(',')
        end = child1[2].index('-',start+1)
        FG_1 = child1[2][start+1:end]
        
        file = open('./molecules_GA/crossover_GA.txt',mode='a',encoding='utf-8')
        file.write('\n'+'Parent1:'+parent1+'\n')
        file.write('Parent2:'+parent2+'\n')
        #file.write(molecule_name+child1[0]+child1[1]+child1[2]+child1[3]+child1[4]+'\n')
        file.close()
        
        #N atom positions
        global position_N_to_add
        position_N_to_add = child1[0]
        position_N_to_add = str(position_N_to_add)
        position_N_to_add = re.findall('[0-9]+',position_N_to_add)
        position_N_to_add = [eval(i) for i in position_N_to_add]
        position_N_to_add = position_N_to_add

        #double bond positions
        global position_db_to_add
        position_db_to_add = child1[1]
        position_db_to_add = str(position_db_to_add)
        position_db_to_add = re.findall('[0-9]+',position_db_to_add)
        position_db_to_add = [eval(i) for i in position_db_to_add]
        position_db_to_add = position_db_to_add

        #fgroup positions
        global position_fg_to_add
        position_fg_to_add = FG_pos_1
        position_fg_to_add = str(position_fg_to_add)
        position_fg_to_add = re.findall('[0-9]+',position_fg_to_add)
        position_fg_to_add = [eval(i) for i in position_fg_to_add]
        position_fg_to_add = position_fg_to_add

        #fgroups
        global fg_to_add
        fg_to_add = FG_1
        fg_to_add = str(fg_to_add)
        fg_to_add = re.findall('[0-9]+',fg_to_add)
        fg_to_add = [eval(i) for i in fg_to_add]
        fg_to_add = fg_to_add

        #linker
        global linker_to_add
        linker_to_add = child1[3]
        linker_to_add = str(linker_to_add)
        linker_to_add = re.findall('[0-9]+',linker_to_add)
        linker_to_add = [eval(i) for i in linker_to_add]
        linker_to_add = linker_to_add

        print('N atoms:', position_N_to_add)
        print('double bonds:', position_db_to_add)
        print('fgroup position:', position_fg_to_add)
        print('fgroup:', fg_to_add)
        print('linker', linker_to_add)
        
        if linker_to_add == [100]:
            double = True
        else:
            double = False
            
        self.double = double
        
        print(self.double)
        
        return position_N_to_add,position_db_to_add,position_fg_to_add,fg_to_add,linker_to_add
    
    def add_n(self):
        
        number_N_to_add = len(position_N_to_add)
        
        self.n_positions = []
        if number_N_to_add == 0:
            position_to_add = ''
            pass
    
        if number_N_to_add > 6:
            number_N_to_add = 6

        if number_N_to_add > 0:
            if self.double == True:
                if 6 in position_N_to_add:
                    position_to_add = position_N_to_add[:-1]
                else:
                    position_to_add = position_N_to_add
            else:
                position_to_add = position_N_to_add

        print('n atoms:',position_to_add)
        
        for i in range(len(position_to_add)):

                #check valency of position to add first, if full, do not add! 
            available = utils_GA.check_valency(position_to_add[i],self.mol_string,self.double)

            self.n_positions.append(position_to_add[i])

            if available == True:
                modified_string = ''
                count_lines_to_position = 0
                for line in self.mol_string.split('\n'):
                    if count_lines_to_position == 3 + position_to_add[i]:
                        line = line.replace(line,line[0:30]+' N  ' + line[34:70])
                    modified_string = modified_string + line + '\n'
                    count_lines_to_position = count_lines_to_position + 1

                self.mol_string = modified_string
    
    def add_fg(self):

        number_fg = len(fg_to_add)

        self.benchmark_atom_number = 6
        self.benchmark_bond_number = 6

        self.fgs = []
        self.save_fg_index = [] 

        for i in range(number_fg):
            
            #the code avoids position 6 if double linker is available, as usually that place has a linker which can be a double bond! 
            if self.double == True:
                if position_fg_to_add == 6:
                    position_to_add = random.choice(1,2,3,4,5)
                else:
                    position_to_add = position_fg_to_add
            else:
                position_to_add = position_fg_to_add
            
            position_to_add = str(position_to_add)[1:-1]
                
            print('fg position:',position_to_add)

                #choose a functional group from  buildingblock files
            random_fg = fg_to_add
            #random_fg = str(random_fg)
            fg_data_file = './buildingblocks_GA/fgroup/%s.mol' %(str(random_fg)[1:-1])
            
            print('fg:',fg_to_add)

                #check availability of position
            available = utils_GA.check_valency(int(position_to_add[i]),self.mol_string,self.double)

                #if available, then do the replacing operation,
            if available == True:

                    #save position and identity of fg, in that order
                self.fgs.append([int(position_to_add[i]),int(random_fg[i])])

                    #build functional group
                self.mol_string, self.benchmark_atom_number, self.benchmark_bond_number, fg_index = utils_GA.build_fg(self.mol_string,fg_data_file,position_to_add,self.benchmark_atom_number,self.benchmark_bond_number)

                self.save_fg_index.append(fg_index)

    def add_db(self):
        
        number_db = len(position_db_to_add)
        
        print('number db:',number_db)

        self.db_positions = []
        
        for i in range(len(position_db_to_add)):
            #double bonds will always be added to forward neighbor (that is if your position 1 then bw 1 and 2, if 6 then between 6 and 1... etc)
            position_for_double = position_db_to_add[i]
            #position_for_double = str(position_for_double)
            #position_for_double = int(position_for_double)

            print('double bonds:',position_for_double)
            
            neighbors = np.zeros((2))
            if position_for_double == 1:
                neighbors[0] = 6
                neighbors[1] = 2
            if position_for_double == 2:
                neighbors[0] = 1
                neighbors[1] = 3
            if position_for_double == 3:
                neighbors[0] = 2
                neighbors[1] = 4
            if position_for_double == 4:
                neighbors[0] = 3
                neighbors[1] = 5
            if position_for_double == 5:
                neighbors[0] = 4
                neighbors[1] = 6
            if position_for_double == 6:
                neighbors[0] = 5
                neighbors[1] = 1

                #3 possible cases: either both neighbors available to double bond, one neighbor, or none

            neighbor1 = int(neighbors[0])
            neighbor2 = int(neighbors[1])

            available_position = utils_GA.check_valency(position_for_double,self.mol_string,self.double)
            available_neighbor = utils_GA.check_valency(neighbor2,self.mol_string,self.double)

            if available_position == True and available_neighbor == True:
                modified_string = ''
                for line in self.mol_string.split('\n'):
                    if ' ' + str(position_for_double) + ' ' in line[0:7] and ' ' + str(neighbor2) + ' ' in line[0:7] and int(line[8]) != 2:
                        line = line.replace(line, line[0:7] + ' 2 ' + line[10:14]) 
                        self.db_positions.append(position_for_double) 

                    modified_string = modified_string + line + '\n'

                self.mol_string = modified_string
    
    def symmetrize(self):

        mol_string = self.mol_string
        number_atoms = self.benchmark_atom_number
        number_bonds = self.benchmark_bond_number
        
        # create a list of atoms that are positions at 2,3,5,6, to switch them,
        # 2 with 6, 3 with 5
        atoms = []
        count_lines = 0
        count_atoms = 1
        for line in self.mol_string.split('\n'):
            if count_lines > 3:
                if count_atoms == 2:
                    if 'N' in line:
                        atoms.append('N')
                    if 'C' in line:
                        atoms.append('C')
                if count_atoms == 6:
                    if 'N' in line:
                        atoms.append('N')
                    if 'C' in line:
                        atoms.append('C')
                if count_atoms == 3:
                    if 'N' in line:
                        atoms.append('N')
                    if 'C' in line:
                        atoms.append('C')
                if count_atoms == 5:
                    if 'N' in line:
                        atoms.append('N')
                    if 'C' in line:
                        atoms.append('C') 
                count_atoms = count_atoms + 1
            
            count_lines = count_lines + 1            

        #switch atom identities 2 and 5, 3 and 5 in string2
        count_lines = 0 
        new_string = ''
        for line in self.mol_string.split('\n'):
            #if you find line 2 switch with the second identity, 6 switch with first
            if count_lines == 3 + 2:
                line = line.replace(atoms[0],atoms[3])
            if count_lines == 3 + 6:
                line = line.replace(atoms[3],atoms[0])
            if count_lines == 3 + 3:
                line = line.replace(atoms[1],atoms[2])
            if count_lines == 3 + 5:
                line = line.replace(atoms[2],atoms[1])
                
            new_string = new_string + line + '\n'
            count_lines = count_lines + 1           
        
        new_string = new_string.strip()
        count_lines = 0
        new_string2 = ''
        for line in new_string.split('\n'):
            if 3+number_atoms < count_lines <= 3+number_atoms + number_bonds and 'END' not in line:
                int_number1 = int(line[1:4])
                int_number2 = int(line[4:7])
                new_number1_int = int_number1 + number_atoms
                new_number2_int = int_number2 + number_atoms
                new_number1 = str(new_number1_int)
                new_number2 = str(new_number2_int)
                if new_number1_int < 10 and new_number2_int < 10:
                    line = line.replace(line[0:7], '  '+ new_number1 + '  ' + new_number2+' ')
                if new_number1_int < 10 and new_number2_int >= 10:
                    line = line.replace(line[0:7], '  ' + new_number1 + ' ' + new_number2+' ')
                if new_number1_int >= 10 and new_number2_int < 10:
                    line = line.replace(line[0:7], ' ' + new_number1 + '  ' + new_number2+' ')
                if new_number1_int >= 10 and new_number2_int >= 10:
                    line = line.replace(line[0:7], ' ' + new_number1 + ' ' + new_number2+' ')
            new_string2 = new_string2 + line + '\n'
            count_lines = count_lines + 1
        
        count_lines = 0
        new_string3 = ''
        for line in new_string2.split('\n'):
            replaced_2 = False
            replaced_3 = False
            if 3+number_atoms < count_lines <= 3+number_atoms + number_bonds:
                if ' ' + str(2+number_atoms) +' ' in line[0:7] or ' ' + str(3+number_atoms) +' ' in line[0:7]:
                    if 2+number_atoms < 10 and 6+number_atoms >= 10:
                        line = line.replace(' ' + str(2+number_atoms) +' ', str(6+number_atoms)+' ')
                        replaced_2 = True
                    else:
                        line = line.replace(' ' + str(2+number_atoms) +' ',' ' + str(6+number_atoms)+' ')
                        replaced_2 = True
                    if 3+number_atoms < 10 and 5+number_atoms >= 10:
                        line = line.replace(' ' + str(3+number_atoms) +' ', str(5+number_atoms)+' ')
                        replaced_3 = True
                    else:
                        line = line.replace(' ' + str(3+number_atoms) +' ',' '+ str(5+number_atoms)+' ')
                        replaced_3 = True
                
                if ' ' + str(5+number_atoms) +' ' in line[0:7] or ' ' + str(6+number_atoms) +' ' in line[0:7]:
                    if replaced_2 == False:
                        if 6+number_atoms >= 10 and 2+number_atoms <10:
                            line = line.replace(str(6+number_atoms) +' ', ' ' + str(2+number_atoms) +' ')
                        else:
                            line = line.replace(' ' + str(6+number_atoms) +' ', ' ' + str(2+number_atoms) +' ')
                    if replaced_3 == False:
                        if 5+number_atoms >=10 and 3+number_atoms < 10:
                            line = line.replace(str(5+number_atoms) +' ', ' ' + str(3+number_atoms) +' ')                        
                        else:
                            line = line.replace(' ' + str(5+number_atoms) +' ', ' ' + str(3+number_atoms) +' ')
            new_string3 = new_string3 + line + '\n'
            count_lines = count_lines + 1

        final_string = ''
        count_lines = 0
        for line in self.mol_string.split('\n'):
            if count_lines <= 3:
                if '  6  6' in line[0:7]:
                    if number_atoms >= 10 and number_bonds >= 10:
                        line = line.replace(' 6  6',str(2*number_atoms)+' '+str(2*number_bonds))
                    if number_atoms >= 10 and number_bonds < 10:
                        line = line.replace(' 6  6',str(2*number_atoms)+' '+str(2*number_bonds))
                    if number_atoms < 10 and number_bonds >= 10:
                        line = line.replace(' 6  6','' + str(2*number_atoms)+' '+str(2*number_bonds))
                    if number_atoms < 10 and number_bonds < 10:
                        line = line.replace(' 6  6','' + str(2*number_atoms)+' '+str(2*number_bonds))
                final_string = final_string + line + '\n'
                
            if 3 < count_lines <= 3 + number_atoms:
                float_x = float(line[2:11])
                new_float_x = float_x - 4.3234233
                new_float_x = round(new_float_x,4)
                new_float_x = str(new_float_x)
                
                line = line.replace(line[3:10],''+new_float_x)
                final_string = final_string + line + '\n'
            count_lines = count_lines + 1
            
            
        count_lines = 0
        for line in new_string3.split('\n'):
            if 3 < count_lines <= 3 + number_atoms:
                final_string = final_string + line + '\n'
            count_lines = count_lines + 1
        
        count_lines = 0
        for line in self.mol_string.split('\n'):
            if  3 + number_atoms < count_lines <= 3+ number_atoms+number_bonds:
                final_string = final_string + line + '\n'
            count_lines = count_lines + 1
            
        count_lines = 0
        for line in new_string3.split('\n'):
            if  3 + number_atoms < count_lines <= 3+ number_atoms+number_bonds:
                final_string = final_string + line + '\n'
            count_lines = count_lines + 1    
        
        final_string = final_string + 'M  END'
        
        self.mol_string = final_string

    def linker(self):

        number_bonds = self.benchmark_bond_number
        number_atoms = self.benchmark_atom_number

        number_bonds = number_bonds*2 + 1
        
        
        if self.double == True:
            random_linker = [100]
            random_linker = str(random_linker)
            linker_data_file = './buildingblocks_GA/linker/%s.mol' %(str(random_linker)[1:-1])

            self.linker = str(random_linker)[1:-1]
        else:
            random_linker = linker_to_add
            random_linker = str(random_linker)
            linker_data_file = './buildingblocks_GA/linker/%s.mol' %(str(random_linker)[1:-1])       

            self.linker = str(random_linker)[1:-1]
            
        #random_linker = str(random_linker)[1:-1]
            
        print('linker:',random_linker)

        count_lines = 0
        modified_string = ''
        for line in self.mol_string.split('\n'):
            if count_lines == 3: 
                line = line.replace(line[0:6],' '+ str(number_atoms*2) + ' ' + str(number_bonds))
            if count_lines == 3+2*number_atoms+1:
                for line2 in fileinput.FileInput(linker_data_file,inplace=0):
                    if 'a' in line2:
                        if number_atoms+6-4 < 10:
                            line2 = line2.replace('a',' '+ str(number_atoms+6-4))
                        else:
                            line2 = line2.replace('a', str(number_atoms+6-4))
                    if 'm' in line2:
                        if number_atoms*2+1 < 10:
                            line2 = line2.replace('m',' ' + str(number_atoms*2+1))
                        elif number_atoms*2+1 >= 10:
                            line2 = line2.replace('m', str(number_atoms*2+1))
                    if 'n' in line2:
                        if number_atoms*2+2 < 10:
                            line2 = line2.replace('n',' ' + str(number_atoms*2+2))
                        elif number_atoms*2+2 >= 10:
                            line2 = line2.replace('n', str(number_atoms*2+2))
                    if 'o' in line2:
                        if number_atoms*2+3 < 10:
                            line2 = line2.replace('o',' ' + str(number_atoms*2+3))
                        elif number_atoms*2+3 >= 10:
                            line2 = line2.replace('o', str(number_atoms*2+3))
                    if 'p' in line2:
                        if number_atoms*2+4 < 10:
                            line2 = line2.replace('p',' ' + str(number_atoms*2+4))
                        elif number_atoms*2+4 >= 10:
                            line2 = line2.replace('p', str(number_atoms*2+4))
                    modified_string = modified_string + line2
                break
            modified_string = modified_string + line + '\n' 
            count_lines = count_lines + 1
        
        modified_string = modified_string + '\n'    
        for line in self.mol_string.split('\n'):
            if ' 1  0 ' in line or ' 2  0 ' in line or ' 3  0 ' in line:
                modified_string = modified_string + line + '\n'
            if 'END' in line:
                modified_string = modified_string + line + '\n'
        
    
        self.mol_string = modified_string


    def finalize(self):
        
        number_atoms = 0
        number_bonds = 0
        count_lines = 0
        for line in self.mol_string.split('\n'):
            if count_lines > 4:
                if 'C' in line or 'N' in line or 'O' in line or 'F' in line or 'Cl' in line or 'S' in line or 'P' in line:
                    number_atoms = number_atoms + 1
                if '1  0' in line[8:12] or '2  0' in line[8:12] or '3  0' in line[8:12]:
                    number_bonds = number_bonds+1
            count_lines = count_lines + 1
        
        final_string = ''
        count_lines = 0
        for line in self.mol_string.split('\n'):
            if '0999' in line:
                line = line.replace(' '+ line[1:3]+' '+line[4:6],' ' + str(number_atoms)+' ' + str(number_bonds))
            final_string = final_string+ line + '\n'
            count_lines = count_lines + 1
            
    #    save_file_path = './final.mol'
        
    #    file = open(save_file_path,mode='w',encoding='utf-8')
        
    #    for line in final_string.split('\n'):
    #        file.write(line+'\n')
    #    file.close()   
        self.number_atoms = number_atoms
        self.number_bonds = number_bonds
        self.mol_string = final_string


    def find_fg(self):

        for file in os.listdir('buildingblocks_GA/fgroup/'):

            fg_lines = []
            for line in fileinput.FileInput('buildingblocks_GA/fgroup/'+file,inplace=0):
                if ' C ' in line or ' N ' in line or ' O ' in line or ' Cl ' in line:
                    fg_lines.append(line[0:71])

            
            mol_string = self.mol_string.split('\n')

            
            fg = True
            for index in range(len(fg_lines)):
                if fg_lines[index] in mol_string:
                    mol_string_index = mol_string.index(fg_lines[index])
                if fg_lines[index] not in mol_string:
                    fg = False
                    break
            
            if fg == True: 
                print('file',file,'index',mol_string_index)


    def remove(self):

        if len(self.fgs) > 0:
            for i in range(len(self.fgs)):
                
                filename = './buildingblocks_GA/fgroup/' + str(self.fgs[i][1]) +'.mol'

                #count number of atoms in functional group from fgroup file in building blockes
                fg_numb_atoms = utils_GA.count_fg_atoms(filename)
                #define the remove indices
                indices_to_remove = [x for x in range(self.save_fg_index[i],self.save_fg_index[i]+fg_numb_atoms)]
                

                new_string = ''
                count_lines = 0 
                shift1 =0
                for line in self.mol_string.split('\n'):
                    
                    if count_lines < 4:
                        new_string = new_string + line + '\n'

                    #first remove the atoms using the indices_to_remove
                    if  3 < count_lines < self.number_atoms + 4 and count_lines-3 not in indices_to_remove:
                        new_string = new_string + line + '\n'  
                    
                    if  self.number_atoms + self.number_bonds + 3 >= count_lines > self.number_atoms + 3 and 'END' not in line:
                        if int(line[1:3]) in indices_to_remove or int(line[4:6]) in indices_to_remove:
                            shift1 = shift1 + 1
                        
                        if int(line[1:3]) not in indices_to_remove and int(line[4:6]) not in indices_to_remove:
                            new_string = new_string + line + '\n'
                    
                    if 'END' in line:
                        new_string = new_string + line + '\n'

                    count_lines = count_lines+1

                #I THINK THE ASSUMPTION THAT THE SHIFT IS ALWAYS THE SAME AS THE NUMBER OF ATOM INDICES TO REMOVE  
                #WILL NOT WORK FOR SOME GROUPS... YOU WILL HAVE TO CALCULATE WHAT THE TRUE BOND SHIFT WOULD BE BY
                #FIRST COUNTING HOW MANY BONDS COME UP USING OUR ATOM INDICES
                

                final_string = utils_GA.adjust(new_string,indices_to_remove[len(indices_to_remove)-1],len(indices_to_remove),self.number_atoms,self.number_bonds, shift1)

                save_mol = open('molecules_GA/removed/mol-files-removed/'+self.molecule_id+'-REM'+str(i)+'.mol',mode='w')
                #save_mol = open('molecules/removed/'+'REM'+str(i)+self.molecule_id+'.mol',mode='w')
                save_mol.write(final_string)
                
                file = open('./molecules_GA/crossover_GA.txt',mode='a',encoding='utf-8')
                file.write('Child1:'+self.molecule_id+'.mol'+'\n')
                file.close()
                

