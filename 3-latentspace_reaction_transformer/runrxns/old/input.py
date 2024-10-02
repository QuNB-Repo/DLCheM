#Dataset
dataset_filepath = 'data/datasets/QM9/qm9.db'


#H-TYPE TRANSFORMATION
targetlabelH_filepath = 'extractembeddings/data/H/model2/qm9-20000/autolab7-2.csv'
#Primary Amine
#H_id = 24
#Secondary Amine
#H_id = 32

#Primary Amide
#H_id = 27

#Primary Alcohol
#H_id = 43
#Secondary Alcohol
#H_id = 40
#Tertiary Alcohol
#H_id = 39

#Carboxylic acids
#H_id = 46

#H-N-=C
#H_id = 35

#H-O-N
#H_id = 1

#new code!!!!!!
#CH3-CH2- Group
#H_id = 17

#OHCH2CH2 --> CH3-O-CH2CH2
#H_id = 69

#1alcseths4order
#H_id = 324

#1alcseths5order
#H_id = 1643
#1alcsethes6order
H_id = 61401

# TARGET EMBEDDING VECTOR TYPE
atom_num = 8

#NOTE number_molecules cannot be in middle of dataset (as the code stands right now)
n_molecules = 20000

#OUTPUT
output_filepath1 = 'fgtransform/data/model2/1alcseths7ordernoopt/init.xyz'
output_filepath2 = 'fgtransform/data/model2/1alcseths7ordernoopt/trans.xyz'
