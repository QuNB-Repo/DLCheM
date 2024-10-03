#!/usr/bin/env python
# coding: utf-8

# In[1]:


import random
import os
import fileinput

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles


# In[2]:


from autocombinatorial_GA import build_GA


# In[3]:


initial_mol_filepath = './buildingblocks_GA/startingmolfile.mol'
max_number_n = 4
max_number_fg = 2
max_number_db = 3
double = True

for i in range(1):
    print('molecule number?')
    molecule_number = input()
    molecule_number = str(molecule_number)
    molecule_name = '%s-' %(molecule_number)
    molecule = build_GA.autocombinatorial_GA(initial_mol_filepath=initial_mol_filepath,max_number_n=max_number_n,max_number_fg=max_number_fg,max_number_db=max_number_db,molecule_name=molecule_name)
    molecule.crossover(double=double)
    molecule.add_n()
    molecule.add_fg()
    molecule.add_db()
    molecule.symmetrize()
    molecule.linker()
    molecule.finalize()
    molecule.show()
    molecule.remove()


# In[4]:


import random
import os
import fileinput

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles

print('molecule number?')
molecule_number = input()
molecule_number = str(molecule_number)

#OPTIMIZING MOL FILE WITH RDKIT
optimized_filename = './molecules_GA/optimized/%s-final-optimized-GA.mol' %(molecule_number)

for file in os.listdir('./molecules_GA/mol-files/mol-files-GA/'):
    if file.startswith('%s-' %(molecule_number)):
        output_file = './molecules_GA/mol-files/mol-files-GA/' + file
        
        mol = Chem.MolFromMolFile(output_file)
        mol = Chem.AddHs(mol)
        molecule=Chem.MolToMolBlock(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol,maxIters=5000)
        Chem.rdmolfiles.MolToMolFile(mol,optimized_filename)

#CONVERTING WITH OBABEL FROM MOL FILE TO GAMESS INPUT FILE
converted_filename = './molecules_GA/converted/%s-final-converted.inp' %(molecule_number)
    
os.system('obabel ' + optimized_filename + ' -O ' + converted_filename)
    
neutgas_filename = './molecules_GA/gamess-files/%s-neutgas-GA.inp' %(molecule_number)
neutdmf_filename = './molecules_GA/gamess-files/%s-neutdmf-GA.inp' %(molecule_number)
dicatdmf_filename = './molecules_GA/gamess-files/%s-dicatdmf-GA.inp' %(molecule_number)
diandmf_filename = './molecules_GA/gamess-files/%s-diandmf-GA.inp' %(molecule_number)

#GAMESS INPUT FOR NEUTRAL GAS PHASE
file = open(neutgas_filename,mode='w',encoding='utf-8')

for line in fileinput.FileInput(converted_filename, inplace=0):
    if '$CONTRL' in line:
        line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=120 ICHARG=0 MULT=1 ISPHER=1 $END' + '\n')
        file.write(line)
        file.write(' $BASIS GBASIS=CCD $END' +'\n')
        file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
        file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
        file.write(' $STATPT HESS=READ $END'+'\n')
for line in fileinput.FileInput(converted_filename, inplace=0):
    if ' $CONTRL COORD=CART' not in line:
        file.write(line)
file.close()

#GAMESS INPUT FOR NEUTRAL DMF SOLVENT
file = open(neutdmf_filename,mode='w',encoding='utf-8')
    
for line in fileinput.FileInput(converted_filename, inplace=0):
    if '$CONTRL' in line:
        line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=120 ICHARG=0 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
        file.write(line)
        file.write(' $BASIS GBASIS=CCD $END' +'\n')
        file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
        file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
        file.write(' $STATPT HESS=READ $END'+'\n')
for line in fileinput.FileInput(converted_filename, inplace=0):
    if ' $CONTRL COORD=CART' not in line:
        file.write(line)
file.close()
    
#GAMESS INPUT FOR DICATION DMF SOLVENT
file = open(dicatdmf_filename,mode='w',encoding='utf-8')
    
for line in fileinput.FileInput(converted_filename, inplace=0):
    if '$CONTRL' in line:
        line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=120 ICHARG=2 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
        file.write(line)
        file.write(' $BASIS GBASIS=CCD $END' +'\n')
        file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
        file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
        file.write(' $STATPT HESS=READ $END'+'\n')
for line in fileinput.FileInput(converted_filename, inplace=0):
    if ' $CONTRL COORD=CART' not in line:
        file.write(line)
file.close()
    
#GAMESS INPUT FOR DIANION DMF SOLVENT
file = open(diandmf_filename,mode='w',encoding='utf-8')
    
for line in fileinput.FileInput(converted_filename, inplace=0):
    if '$CONTRL' in line:
        line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=120 ICHARG=-2 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
        file.write(line)
        file.write(' $BASIS GBASIS=CCD $END' +'\n')
        file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
        file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
        file.write(' $STATPT HESS=READ $END'+'\n')
for line in fileinput.FileInput(converted_filename, inplace=0):
    if ' $CONTRL COORD=CART' not in line:
        file.write(line)
file.close()
        
#os.system('rm *mol')
#os.system('rm *opt*')
os.system('rm *converted*')


# In[5]:


for file_a in os.listdir('./molecules_GA/removed/mol-files-removed'):
    if file_a.endswith('-REM0.mol'):
        if file_a.startswith('%s-' %(molecule_number)):
            output_a_file = './molecules_GA/removed/mol-files-removed/' + file_a
            print(output_a_file)
            optimized_a_filename = './molecules_GA/removed/optimized-removed/%s-a-removed-optimized.mol' %(molecule_number)

            mol = Chem.MolFromMolFile(output_a_file)
            mol = Chem.AddHs(mol)
            molecule=Chem.MolToMolBlock(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol,maxIters=5000)
            Chem.rdmolfiles.MolToMolFile(mol,optimized_a_filename)
            
            converted_a_filename = './molecules_GA/removed/converted-removed/%s-a-removed-converted.inp' %(molecule_number)

            os.system('obabel ' + optimized_a_filename + ' -O ' + converted_a_filename)

            neutdmf_a_filename = './molecules_GA/removed/gamess-files-removed/%s-a-neutdmf-GA.inp' %(molecule_number)
            dicatdmf_a_filename = './molecules_GA/removed/gamess-files-removed/%s-a-dicatdmf-GA.inp' %(molecule_number)
            diandmf_a_filename = './molecules_GA/removed/gamess-files-removed/%s-a-diandmf-GA.inp' %(molecule_number)

            #GAMESS INPUT FOR NEUTRAL DMF SOLVENT
            file = open(neutdmf_a_filename,mode='w',encoding='utf-8')

            for line in fileinput.FileInput(converted_a_filename, inplace=0):
                if '$CONTRL' in line:
                    line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=120 ICHARG=0 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
                    file.write(line)
                    file.write(' $BASIS GBASIS=CCD $END' +'\n')
                    file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
                    file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
                    file.write(' $STATPT HESS=READ $END'+'\n')
            for line in fileinput.FileInput(converted_a_filename, inplace=0):
                if ' $CONTRL COORD=CART' not in line:
                    file.write(line)
            file.close()

            #GAMESS INPUT FOR DICATION DMF SOLVENT
            file = open(dicatdmf_a_filename,mode='w',encoding='utf-8')

            for line in fileinput.FileInput(converted_a_filename, inplace=0):
                if '$CONTRL' in line:
                    line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=120 ICHARG=2 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
                    file.write(line)
                    file.write(' $BASIS GBASIS=CCD $END' +'\n')
                    file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
                    file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
                    file.write(' $STATPT HESS=READ $END'+'\n')
            for line in fileinput.FileInput(converted_a_filename, inplace=0):
                if ' $CONTRL COORD=CART' not in line:
                    file.write(line)
            file.close()

            #GAMESS INPUT FOR DIANION DMF SOLVENT
            file = open(diandmf_a_filename,mode='w',encoding='utf-8')

            for line in fileinput.FileInput(converted_a_filename, inplace=0):
                if '$CONTRL' in line:
                    line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=120 ICHARG=-2 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
                    file.write(line)
                    file.write(' $BASIS GBASIS=CCD $END' +'\n')
                    file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
                    file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
                    file.write(' $STATPT HESS=READ $END'+'\n')
            for line in fileinput.FileInput(converted_a_filename, inplace=0):
                if ' $CONTRL COORD=CART' not in line:
                    file.write(line)
            file.close()


# ###### 
