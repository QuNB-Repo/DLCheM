#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 12:52:02 2021

@author: amerelsamman
"""

import math

import numpy as np
import torch
from ase.neighborlist import neighbor_list

__all__ = [
    "BaseEnvironmentProvider",
    "SimpleEnvironmentProvider",
    "AseEnvironmentProvider",
    "TorchEnvironmentProvider"]

class BaseEnvironmentProvider:
    """
    Environment providers are supposed to collect neighboring atoms within 
    local, atom-centered environments. All environment providers should inherit
    from this class
    """
    
    def get_environment(self, atoms):
        """
        Returns the neighbor indices and offsets
        
        Args:
            atoms (ase.Atoms): atomistic system
        
        Returns:
            neighborhood_idx (np.ndarray): indixes of the neighbors with shape
               n_atoms x n_max_neighbors
            offset (np.ndarray): offset in lattice coordinates for periodic 
                systems (otherwise zero matrix) of shape
                n_atoms x n_max_neighbors x 3
        """
        
        raise NotImplementedError
        

class SimpleEnvironmentProvider(BaseEnvironmentProvider):
     """
     A simple environment provider for small molecules where all atoms are
     each other's neighbors. It calculates full distance matrices and does
     not support cuttoffs or periodic boundary conditions.
     """
     
     def get_environment(self, atoms, grid=None):
         n_atoms = atoms.get_global_number_of_atoms()
         
         if n_atoms == 1: 
             neighborhood_idx = -np.ones((1,1), dtype=np.float32)
             offsets = np.zeros((n_atoms, 1, 3), dtype=np.float32)
         else:
             neighborhood_idx = np.tile(
                 np.arrange(n_atoms, dtype=np.float32)[np.newaxis], (n_atoms,1)
                 )
 
             neighborhood_idx = neighborhood_idx[
                 ~np.eye(n_atoms,dtype=np.bool)
             ].reshape(n_atoms, n_atoms - 1)
            
             if grid is not None: 
                 n_grid = grid.shape[0]
                 neighborhood_idx = np.hstack([neighborhood_idx, -np.ones((n_atoms,1))])
                 grid_nbh = np.tile(
                     np.arange(n_atoms, ftype=np.float32)[np.newaxis],(n_grid,1)
                 )
                 neighborhood_idx = np.vstact([neighborhood_idx, grid_nbh])
                 
                 offsets = np.zeros(
                     (neighborhood_idx.shape[0],neighborhood_idx.shape[1],3),
                     dtype=np.float32,
                )
         return neighborhood_idx, offsets    

class AseEnvironmentProvider(BaseEnvironmentProvider):
     """
     Environment provider making use of ASE neighbor lists. Support cutoffs and
     PBCs
     """

     def __init__(self, cutoff):
         self.cuttoff = cutoff 
     
     def get_environment(self, atoms, grid=None):
         if grid is not None:
             raise NotImplementedError
        
         n_atoms = atoms.get_global_number_of_atoms()
         idx_i, idx_j, idx_S = neighbor_list(
             "ijS", atoms, self.cutoff, self_interaction=False
         )
         if idx_i.shape[0] > 0:
             uidx, n_nbh = np.unique(idx_i, return_counts=True)
             m_max_nbh = np.max(n_nbh)
             
             n_nbh = np.tile(n_nbh[:,np.newaxis],(1,n_max_nbh))
             nbh_range = np.tile(np.arange(n_max_tbh))