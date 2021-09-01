# -*- coding: utf-8 -*-
"""

"""

from ase.io import read
import numpy as np
from schnetpack import AtomsData


def generate(datasets_file_dir,available_properties,index):

    index = str(index)
    atoms = read(datasets_file_dir  + '.xyz', index=':'+index)
    
    property_list = []
    
    for at in atoms:
        energy = np.array([float(list(at.info.keys())[0])], dtype=np.float32)
        property_list.append(
            {'energy': energy}
        )
        
    print('Properties:', property_list)
    
    new_dataset = AtomsData(datasets_file_dir  + '.db', available_properties=[available_properties])
    new_dataset.add_systems(atoms,property_list)
    
