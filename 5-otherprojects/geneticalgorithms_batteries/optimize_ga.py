import os
from rdkit import Chem
from rdkit.Chem import AllChem
import fileinput


def optimize_GA_mol(MOLECULE_NUMBER,MOLECULES_FILEPATH,OUTPUT_FILEPATH,INP_CONVERTED_FILEPATH,GAMESS_FILEPATH):
    '''
    Writing the GAMESS-US file for genetic algorithm constructed molecules so that 
    We can carry out computation on its using GAMESS-US and the specified level
    of electronic structure theory B3LYP/CCD

        MOLECULE_NUMBER 
            - Molecule number of the GA file to optimize and create GAMESS
            input file 
        MOLECULES_FILEPATH
            - filepath where the genetically produced molecule files are stored
        OUTPUT_FILEPATH
            - filepath where the optimize GA molecule file is to be stored for all types
            neutral-gas, neutral-dmf, cation-gas, cation-dmf, anion-gas, anion-dmf
        INP_CONVERTED_FILEPATH
            - filepath to store .inp converted molecules (using obabel) 
            to be used for quick visualization
        GAMESs_FILEPATH
            - filepath to store the final GAMESS input files for all oxidation states of bpy derivative:
            cation, anion, neutral, and in both gas and DMF solvent. Level of theory is B3LYP/CCD
        
    
    '''

    print('molecule number?')
    molecule_number = str(MOLECULE_NUMBER)

    #OPTIMIZING MOL FILE WITH RDKIT
    optimized_filename = OUTPUT_FILEPATH + '%s-final-optimized-GA.mol' %(molecule_number)

    for file in os.listdir(MOLECULES_FILEPATH):
        if file.startswith('%s-' %(molecule_number)):
            mol_filepath = MOLECULES_FILEPATH + file
            
            mol = Chem.MolFromMolFile(mol_filepath)
            mol = Chem.AddHs(mol)
            molecule=Chem.MolToMolBlock(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol,maxIters=5000)
            Chem.rdmolfiles.MolToMolFile(mol,optimized_filename)

    #CONVERTING WITH OBABEL FROM MOL FILE TO GAMESS INPUT FILE
    converted_filename = INP_CONVERTED_FILEPATH + '%s-final-converted.inp' %(MOLECULE_NUMBER)
        
    os.system('obabel ' + optimized_filename + ' -O ' + converted_filename)
        
    neutgas_filename = GAMESS_FILEPATH + '%s-neutgas-GA.inp' %(molecule_number)
    neutdmf_filename = GAMESS_FILEPATH + '/%s-neutdmf-GA.inp' %(molecule_number)
    dicatdmf_filename = GAMESS_FILEPATH + '%s-dicatdmf-GA.inp' %(molecule_number)
    diandmf_filename = GAMESS_FILEPATH + '%s-diandmf-GA.inp' %(molecule_number)

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
