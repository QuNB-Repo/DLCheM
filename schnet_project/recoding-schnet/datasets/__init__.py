#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 12:17:47 2021

@author: amerelsamman
"""


from schnetpack import AtomsData, get_center_of_mass
from schnetpack.data.atoms import logger 
from schnetpack.environment import SimpleEnvironmentProvider

class DownloadableAtomsData(AtomsData):
    '''
    Base class for online datasets that can be automatically downloaded
    Args:
        dbpath (str): path to directory containing database