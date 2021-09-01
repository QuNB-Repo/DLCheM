# -*- coding: utf-8 -*-
"""

"""

import fileinput
from schnetpack.datasets import QM9

from schnetpack import AtomsData

dataset_file_path = '../../../data/datasets/extendchain/extendchain.db'
data = AtomsData(dataset_file_path,available_properties=['energy'])
idx=1
gamess_file_path = '../../../data/qm9_to_gamess/%s.inp' %(idx)
template_path = '../../../data/qm9_to_gamess/template.inp'
basis_file_dir = '../../../data/qm9_to_gamess/basis/631g2dfpi/'


at, props = data.get_properties(idx)

number_atoms = len(props['_atomic_numbers'])


lines = []
for i in range(number_atoms):
    x = props['_positions'][i][0]
    y = props['_positions'][i][1] 
    z = props['_positions'][i][2]
    x = x.detach().numpy()
    y = y.detach().numpy()
    z = z.detach().numpy()
    if props['_atomic_numbers'][i] == 1:
        lines.append([])
        element_line = 'H 1.0 ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
        basis_file = fileinput.FileInput(basis_file_dir+'H.inp',inplace=0)
        basis_lines = ([])
        count=0
        for line in basis_file:
            basis_lines.append([])
            basis_lines[count] = line
            count=count+1
        lines[i] = element_line
        for j in range(count):
            lines[i] = lines[i] + basis_lines[j] 
    if props['_atomic_numbers'][i] == 6:
        lines.append([])
        element_line = 'C 6.0 ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
        basis_file = fileinput.FileInput(basis_file_dir+'C.inp',inplace=0)
        basis_lines = ([])
        count=0
        for line in basis_file:
            basis_lines.append([])
            basis_lines[count] = line
            count=count+1
        lines[i] = element_line
        for j in range(count):
            lines[i] = lines[i] + basis_lines[j] 
    if props['_atomic_numbers'][i] == 7:
        lines.append([])
        element_line = 'N 7.0 ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
        basis_file = fileinput.FileInput(basis_file_dir+'N.inp',inplace=0)
        basis_lines = ([])
        count=0
        for line in basis_file:
            basis_lines.append([])
            basis_lines[count] = line
            count=count+1
        lines[i] = element_line
        for j in range(count):
            lines[i] = lines[i] + basis_lines[j] 
    if props['_atomic_numbers'][i] == 8:
        lines.append([])
        element_line = 'O 8.0 ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
        basis_file = fileinput.FileInput(basis_file_dir+'O.inp',inplace=0)
        basis_lines = ([])
        count=0
        for line in basis_file:
            basis_lines.append([])
            basis_lines[count] = line
            count=count+1
        lines[i] = element_line
        for j in range(count):
            lines[i] = lines[i] + basis_lines[j] 
    if props['_atomic_numbers'][i] == 9:
        lines.append([])
        element_line = 'F 9.0 ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
        basis_file = fileinput.FileInput(basis_file_dir+'F.inp',inplace=0)
        basis_lines = ([])
        count=0
        for line in basis_file:
            basis_lines.append([])
            basis_lines[count] = line
            count=count+1
        lines[i] = element_line
        for j in range(count):
            lines[i] = lines[i] + basis_lines[j] 

file = open(gamess_file_path,mode='w',encoding='utf-8')

template = fileinput.FileInput(template_path,inplace=0)
for line in template:
   file.write(line) 

for i in range(len(lines)):
    string = lines[i]
    file.write(string) 
    
file.write(' $END' + '\n')
file.close()