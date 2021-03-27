#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 06:53:33 2020

@author: amerelsamman
"""

import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import LSTM
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import Dense


#construct dictionary of SMILES


symbols =	{
  "/": 0,
  "=": 1,
  "C": 2,
  "H": 3,
}
print(symbols)

data = ('HCH')
print(data)

model = keras.Sequential()
model.add(layers.LSTM(128))
model.add(LSTM(128, input_shape=(data.shape[1], data.shape[2]), return_sequences = True))
model.add(Dropout(0.2))
model.add(LSTM(256, return_sequences = True))
model.add(Dropout(0.2))
model.add(LSTM(512, return_sequences = True))
model.add(Dropout(0.2))
model.add(LSTM(256, return_sequences = True))
model.add(Dropout(0.2))
model.add(LSTM(128))
model.add(Dropout(0.2))
model.add(Dense(Y.shape[1], activation='softmax'))