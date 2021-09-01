# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 14:05:20 2021

@author: aelsamma
"""

import pandas as pd

import networkx as nx




data = pd.read_csv('../../../data/network_test3_7.csv', delimiter=',')


G = nx.from_pandas_edgelist(data,
                            source = 'Source',
                            target = 'Target',
                            edge_attr = 'Weight')

nx.draw_networkx(G)

