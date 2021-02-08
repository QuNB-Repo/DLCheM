#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 11:56:33 2021

@author: amerelsamman
"""


import logging
import os
import re
import shutil
import tarfile
import tempfile
from urllib import request as request

import numpy as np 
from ase.io.extxyz import read_xyz
from ase.units import Debye, Bohr, Hartree, eV

import schnetpack as spk
from schnetpack.datasets import DownloadableAtomsData

