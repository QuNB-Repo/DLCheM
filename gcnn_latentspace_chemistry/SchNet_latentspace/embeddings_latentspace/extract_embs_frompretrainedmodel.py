from schnetpack.datasets import QM9
import schnetpack as spk
import torch 

import numpy as np
from numpy import genfromtxt, savetxt

from schnetpack import AtomsData

from label.manuallabel2 import labeller
from label.manuallabel2.utils import utils


'''
Last Updated: 2024-03-08
'''

def extract_embs(QM9_TRUE,DB_FILEPATH, MODEL_FILEPATH,SAVE_FILEPATH,START,END,N_FEATURES,LAYERS,ELEMENTS=[1,6,7,8,9],AVAILABLE_PROPERTIES='N/A',LABEL=True,ADD_HEADER=True,RESTRICT_LABEL=False,ALLOWED_LABELS=[],INDICES_TO_LABEL=[]):
    '''
    
    Extracting & labeling embeddings of QM9 using a pretrained SchNet model
    and indices list for each molecule (a list of lists) that point to which atoms to extract

    Args:
        QM9_TRUE                - a boolean depending on whether QM9.db is being used because we
                                   load all the properties using QM9 method instead
        DB_FILEPATH             - the file location for the db file to exctract embeddings of
        MODEL_FILEPATH          - the file location for trained model (best_model)
        SAVE_FILEPATH           - the file location saving the embeddings (and interaction/convolution layers)
        START                   - the start index of the db file to extract embs
        END                     - the END index of the db file to extract embs
        N_FEATURES              - the number of features of the embedding in the trained model
        LAYERS                  - the layer range to extract embeddings
        ELEMENTS
        AVAILABLE_PROPERTIES    - if loading a different dataset than QM9, need to provide its available
                                   properties as a list of strings ['property1','property2',...]  
        LABEL                   - boolean to decide whether to LABEL the extracted embedding with manual LABEL code
        ADD_HEADER              - boolean whether to add a header to the resulting embedding csv files
        RESTRICT_LABEL
        ALLOWED_LABELS          -
        INDICES_TO_LABEL        - A lists of lists, one list for each molecule, informs algorithm which atom indices
                                  to extract
    
    Processing:
        
        db_data                 - the db data (molecule xyz and properties), from DATASET_FILEPATH
        device                  - device used to access model (can use either CPU or CUDA)
        model                   - loading the pretrained schnet model using torch, from MODEL_FILEPATH
        converter               - setting the atomsconverter which converts at to inputs ready for schnet
                                  USES atomic simulation environment (ase) to extract data about the molecule (defined in at)
        save_filepathvs         - defined as the same as SAVE_FILEPATH but replacing '.csv' with 'vs.csv' 
                                  it is the associated interaction updates for the embeddings
        save_filepathws         - defined as the same as SAVE_FILEPATH but replacing '.csv' with 'ws.csv'
                                  contains the filter weights for each atom 
                                  (the size of which depends on number of neighbors in molecule and N_FEATURES size)   
        xs_header               - created header for embeddings 
                                  'embs1,embs2,embs3,...,embs_N_FEATURES,mol_index,atom_index,atomic_number,coordx,coordy,coordz,ldalabel,gnumarker,gnucolor,hexcolor,layer,atomizationE,fg_key \n'
        vs_header               - created header for interactions updates
                                   'vs1,vs2,vs3,...,vs_N_FEATURES,mol_index,atom_index,atomic_number,coordx,coordy,coordz,ldalabel,gnumarker,gnucolor,hexcolor,layer,atomizationE,fg_key \n'   
    
    Returns:
        save_xs_datafile        - write the extracted embeddings for layers specified 
                                   of each atom and labels (if needed) into a file
        save_vs_datafile        - write the interaction updates for layers specified
                                   of each atom and labels (if needed) 
        save_ws_datafile        - write the weight filters for layers specified 
                                   of each atom and labels (if needed) into a file 
    '''

    #Load testing data, if QM9 use the qm9 loader, otherwise, use atomsdata loader and use available properties
    if QM9_TRUE == True:
        db_data = QM9(DB_FILEPATH,download=False,remove_uncharacterized=True)
    else:
        db_data = AtomsData(DB_FILEPATH,available_properties=AVAILABLE_PROPERTIES)

    #load model
    device = 'cpu'
    model = torch.load(MODEL_FILEPATH,map_location=torch.device(device))

    # load the data converter which will convert the data to machine-readable format for the algorithm
    converter = spk.data.AtomsConverter(device=device)

    #open a save file
    save_xs_datafile = open(SAVE_FILEPATH,mode='a')
    save_filepathvs = SAVE_FILEPATH.replace('.csv','v.csv')
    save_vs_datafile = open(save_filepathvs,mode='a')
    save_filepathWs = SAVE_FILEPATH.replace('.csv','W.csv')
    save_Ws_datafile = open(save_filepathWs,mode='a')

    #write header of the csv file
    if ADD_HEADER == True:
        xs_header = ','.join('embs'+str(x) for x in range(N_FEATURES)) +','+'mol_index,atomic_number,coordx,coordy,coordz,ldalabel,gnumarker,gnucolor,hexcolor,layer,atomizationE,fg_key \n'
        vs_header = ','.join('vs'+str(x) for x in range(N_FEATURES)) +','+'mol_index,atomic_number,coordx,coordy,coordz,ldalabel,gnumarker,gnucolor,hexcolor,layer,atomizationE,fg_key \n'
        save_xs_datafile.write(xs_header)
        save_vs_datafile.write(vs_header) 
    
    #for each_molecule in the start --> end of db dataset
    for each_molecule in range(START,END):


        #tiny load bar
        if each_molecule % 1000 == 0:
            print(each_molecule)

        print(each_molecule)
        #get properties out of the db_data for each molecule
        at, props = db_data.get_properties(each_molecule)

        #convert to machine readable format
        inputs = converter(at)

        #use model on inputs  
        model(inputs)

        #define positions and atomic numbers and number_atoms
        positions = props['_positions'].detach().numpy()
        atomic_numbers = props['_atomic_numbers'].detach().numpy()
        number_atoms = len(atomic_numbers)

        if INDICES_TO_LABEL == []:
            indices_to_label = [each_mol_idx for each_mol_idx in range(number_atoms)]
        else:
            indices_to_label = [INDICES_TO_LABEL[each_molecule]]
    

        from schnetpack.representation.schnet import xs, vs, Ws
        from schnetpack.atomistic.output_modules import yi

        for layer in range(LAYERS[0],LAYERS[1]):
            xemb = xs[layer+1].cpu().detach().numpy()
            vsadd = vs[layer].cpu().detach().numpy()
            Wslayer = Ws[layer].cpu().detach().numpy()

            #write xyz and mol file in a temp file
            mol_filename = utils.xyz2mol(props)

            #for each atom in the molecule
            for each_atom in range(number_atoms):
                if each_atom in indices_to_label:
                    if atomic_numbers[each_atom] in ELEMENTS:
                        if LABEL == True:
                            #run labeller code to output LDA label, gnuplot color label, gnuplot marker label... NEXT UPDATE!!! 
                            fg_key,ldalabel,gnucolor,gnumarker,hexcolor= labeller.label(mol_filename,number_atoms,each_atom,atomic_numbers[each_atom],each_molecule)
                        else:
                            fg_key = ''
                            ldalabel = 0
                            gnucolor = 0
                            gnumarker = 0
                            hexcolor = ''

                        #write a row of information for the csv file, the row will contaion
                        #molecule_number, the atomic number of atom, ae prediciton, n_features embedding[layer],
                        #same for the n_features of vs[layer],
                        if RESTRICT_LABEL == True:
                            fg_key,ldalabel,gnucolor,gnumarker,hexcolor= labeller.label(mol_filename,number_atoms,each_atom,atomic_numbers[each_atom],each_molecule)
                            if ldalabel in ALLOWED_LABELS:
                                row_embs = ','.join(str(x) for x in xemb[0][each_atom][:N_FEATURES]) +','+str(each_molecule) + ','+ str(atomic_numbers[each_atom]) +',' + ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+ ','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
                                row_vs = ','.join(str(x) for x in vsadd[0][each_atom][:N_FEATURES]) +','+str(each_molecule) + ','+ str(atomic_numbers[each_atom]) +','+ ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+ ','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
                                save_xs_datafile.write(row_embs)
                                save_vs_datafile.write(row_vs)
                        else:
                            row_embs  = ','.join(str(x) for x in xemb[0][each_atom][:N_FEATURES]) +','+str(each_molecule) + ','+ str(atomic_numbers[each_atom]) +',' + ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
                            row_vs    = ','.join(str(x) for x in vsadd[0][each_atom][:N_FEATURES]) +','+str(each_molecule) + ','+ str(atomic_numbers[each_atom]) +','+ ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+ ','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
#                            flatmatrix_Ws = ','.join(str(x) for x in Wslayer[0][each_atom][:number_atoms-1][:n_features]) + ',' + str(atomic_numbers[each_atom]) + ',' + ','.join(str(x) for x in positions[each_atom][:3])+','+str(ldalabel)+','+str(gnumarker)+','+str(gnucolor)+','+str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom].detach().numpy()) + ',' + fg_key + '\n'

                            row_embs = ','.join(str(x) for x in xemb[0][each_atom][:N_FEATURES]) +','+str(each_molecule) + ','+ str(atomic_numbers[each_atom]) +',' + ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+ ','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
                            row_vs = ','.join(str(x) for x in vsadd[0][each_atom][:N_FEATURES]) +','+str(each_molecule) + ','+ str(atomic_numbers[each_atom]) +','+ ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+ ','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
                            save_xs_datafile.write(row_embs)
                            save_vs_datafile.write(row_vs)
#                            save_Ws_datafile.write(flatmatrix_Ws)


    save_xs_datafile.close()
    save_vs_datafile.close()


