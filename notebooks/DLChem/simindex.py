# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 19:00:06 2021

@author: aelsamma
"""

import pandas as pd
import os

import rdkit 

from schnetpack.datasets import QM9

from utils import utils


#label file
element = 'O'
name_label = '5000'
label_file = '../../../data/label_qm9/O/label%s%s.csv' %(element,name_label)

#data file
qm9data = QM9('../../../data/QM9/qm9.db', download=False, remove_uncharacterized=True)

number_molecules = 1

for idx in range(number_molecules):
    
    at, props = qm9data.get_properties(idx)
    
    utils.write_xyz_from_db(props,idx)
    
    utils.xyz_to_mol(idx)

    

label = pd.read_csv(label_file,delimiter=',')




