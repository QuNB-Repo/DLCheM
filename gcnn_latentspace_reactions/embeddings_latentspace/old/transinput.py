element = 'O'
atom_num = 8
db_file_path = 'fgtransform/data/model2/test/init.db'

use_label = False
#label_filepath = './fgtransform/data/model2/2alcseths/init161Olabel.csv'
#label_filepath = './trans32Nlabel.csv'
#label_id = 4

#number of molecules in the db
number_inputs=2

#to load checkpoint
qm9_file = 'data/datasets/QM9/qm9.db'

save_filepath = 'fgtransform/data/model2/test/init2Oembs.csv'
newlabel = 0 

split_file='data/trainedmodels/model2/split.npz' 
checkpoint_path = 'data/trainedmodels/model2/trained.pth' 
n_atom_basis=30
n_filters=30
n_gaussians=20
n_interactions=3
cutoff = 4.
index = number_inputs
