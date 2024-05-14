
dataset_filepath = 'data/datasets/QM9/qm9.db'

targetlabelH_filepath = 'extractembeddings/data/H/model2/qm9-20000/Hlabel0-20000.csv'

#primary alcohol H's for oxidation
H1_id = 43
H2_id = 6

#secondary
#H1_id = 40
#H2_id = 10

#target embedding vector type
atom_num = 8

#number of molecules, CANNOT BE IN MIDDLE OF DATASET AS CODE STANDS. CODE REQUIRES COUNTING ATOMS PER MOLECULE
n_molecules = 20000

output_filepath1 = 'fgtransform/data/model2/1alcsalds/init.xyz'
output_filepath2 = 'fgtransform/data/model2/1alcsalds/trans.xyz'

