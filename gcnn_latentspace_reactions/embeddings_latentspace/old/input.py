element = 'N'
atom_num = 7
db_file_path = 'fgtransform/data/model2/1amns2amns/trans.db'

use_label = True
#label_filepath = './trans83Olabel.csv'
label_filepath = './trans32Nlabel.csv'
label_id = 12

#number of molecules in the db
number_inputs=32

#to load checkpoint
qm9_file = 'data/datasets/QM9/qm9.db'

save_filepath = 'trans32Nembs.csv'

split_file='data/trainedmodels/model2/split.npz' 
checkpoint_path = 'data/trainedmodels/model2/trained.pth' 
n_atom_basis=30
n_filters=30
n_gaussians=20
n_interactions=3
cutoff = 4.
index = number_inputs
