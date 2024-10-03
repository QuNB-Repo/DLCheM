
#initial code which looks through all the mol files will just look for the symbol 'N-' 
import os


class remove_fg():
    def __init__(self,mol_dir=''):
        
        filenames = os.listdir(mol_dir)

        mol_files = []
        for file in filenames:
            if 'N' in file and 'DB' in file and '.mol' in file:
                mol_files.append(file)

        self.mol_files = mol_files
    

    def remove()