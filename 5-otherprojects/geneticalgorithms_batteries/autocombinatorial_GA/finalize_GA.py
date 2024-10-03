import random
import os
import fileinput

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles

def optimize(self,molecule_number,output_filename,optimized_filename):
    
    output_filename = './molecules_GA/%s-' + self.molecule_id + '.mol' %(molecule_number)
    optimized_filename = './molecules_GA/optimized/%s-final-optimized.mol' %(molecule_number)
    
    mol = Chem.MolFromMolFile(output_filename)
    mol = Chem.AddHs(mol)
    molecule=Chem.MolToMolBlock(mol)
    AllChem.MMFFOptimizeMolecule(mol,maxIters=5000)
    Chem.rdmolfiles.MolToMolFile(mol,optimized_filename)
    
def convert(self,molecule_number,optimized_filename,converted_filename,trim_filename,neutgas_filename,neutdmf_filename,dicatdmf_filename,diandmf_filename):
    
    converted_filename = './molecules_GA/converted/%s-final-converted.inp' %(molecule_number)
    trim_filename = './molecules_GA/%s-final-trim.inp' %(molecule_number)
    neutgas_filename = './molecules_GA/neutral-gas/%s-neutgas.inp' %(molecule_number)
    neutdmf_filename = './molecules_GA/neutral-dmf/%s-neutdmf.inp' %(molecule_number)
    dicatdmf_filename = './molecules_GA/dication-dmf/%s-dicatdmf.inp' %(molecule_number)
    diandmf_filename = './molecules_GA/dianion-dmf/%s-diandmf.inp' %(molecule_number)
    
    os.system('obabel ' + optimized_filename + ' -O ' + converted_filename)

#GAMESS INPUT FOR NEUTRAL GAS PHASE
    file = open(neutgas_filename,mode='w',encoding='utf-8')

    for line in fileinput.FileInput(converted_filename, inplace=0):
        if '$CONTRL' in line:
            line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=30 ICHARG=0 MULT=1 ISPHER=1 $END' + '\n')
            file.write(line)
   # final_filename.write(line)
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
            line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=30 ICHARG=0 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
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
            line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=80 ICHARG=2 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
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
            line = line.replace(line, ' $CONTRL SCFTYP=RHF RUNTYP=HESSIAN EXETYP=RUN DFTTYP=B3LYP' + '\n' 'MAXIT=80 ICHARG=-2 MULT=1 ISPHER=1 $END' + '\n' ' $PCM SOLVNT=DMF SMD=.TRUE. $END' + '\n' ' $FORCE METHOD=SEMINUM $END' + '\n')
            file.write(line)
            file.write(' $BASIS GBASIS=CCD $END' +'\n')
            file.write(' $SYSTEM TIMLIM=1300 MWORDS=500 PARALL=.TRUE. $END'+'\n')
            file.write(' $SCF DIRSCF=.TRUE. $END '+'\n')
            file.write(' $STATPT HESS=READ $END'+'\n')
    for line in fileinput.FileInput(converted_filename, inplace=0):
        if ' $CONTRL COORD=CART' not in line:
            file.write(line)
    file.close()