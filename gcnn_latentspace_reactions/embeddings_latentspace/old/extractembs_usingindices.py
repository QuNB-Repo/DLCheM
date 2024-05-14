from schnetpack.datasets import QM9
import schnetpack as spk
import torch 

import numpy as np
from numpy import genfromtxt, savetxt

from schnetpack import AtomsData

from LABEL.manuallabel2 import labeller
from LABEL.manuallabel2.utils import utils

'''
Last Updated: 2024-03-01
'''

def extractembs(QM9,DATASET_FILEPATH,MODEL_FILEPATH,SAVE_FILEPATH,start,END,N_FEATURES,LAYER_RANGE,AVAILABLE_PROPERTIES,LABEL,ADD_HEADER,INDICES_TO_LABEL):
    '''
    
    Extracting & labeling embeddings of QM9 using a pretrained SchNet model
    and indices list for each molecule that point to which atoms to extract

    Args:
        QM9                     - a boolean depending on whether QM9.db is being used because we
                                   load all the properties using QM9 method instead
        DATASET_FILEPATH        - the file location for the db file to exctract embeddings of
        MODEL_FILEPATH          - the file location for trained model (best_model)
        SAVE_FILEPATH           - the file location saving the embeddings (and interaction/convolution layers)
        start                   - the start index of the db file to extract embs
        END                     - the END index of the db file to extract embs
        N_FEATURES              - the number of features of the embedding in the trained model
        LAYER_RANGE             - the layer range to extract embeddings
        AVAILABLE_PROPERTIES    - if loading a different dataset than QM9, need to provide its available
                                   properties as a list of strings ['property1','property2',...]  
        LABEL                   - boolean to decide whether to LABEL the extracted embedding with manual LABEL code
        ADD_HEADER              - boolean whether to add a header to the resulting embedding csv files
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
    #Load testing data, if QM9 use the QM9 loader, otherwise, use atomsdata loader and use available properties
    if QM9 == True:
        db_data = QM9(DATASET_FILEPATH,download=False,remove_uncharacterized=True)
    else:
        db_data = AtomsData(DATASET_FILEPATH,AVAILABLE_PROPERTIES=AVAILABLE_PROPERTIES)

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
    save_ws_datafile = open(save_filepathWs,mode='a')

    #write header of the csv file
    if ADD_HEADER == True:
        xs_header = ','.join('embs'+str(x) for x in range(N_FEATURES)) +','+'mol_index,atom_index,atomic_number,coordx,coordy,coordz,ldalabel,gnumarker,gnucolor,hexcolor,layer,atomizationE,fg_key \n'
        vs_header = ','.join('vs'+str(x) for x in range(N_FEATURES)) +','+'mol_index,atom_index,atomic_number,coordx,coordy,coordz,ldalabel,gnumarker,gnucolor,hexcolor,layer,atomizationE,fg_key \n'
        save_xs_datafile.write(xs_header)
        save_vs_datafile.write(vs_header) 
    
    #for each_molecule in the start --> END of db dataset
    for each_molecule in range(start,END):


        #tiny load bar
        if each_molecule % 1000 == 0:
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

        from schnetpack.representation.schnet import xs, vs, Ws
        from schnetpack.atomistic.output_modules import yi

        for layer in range(LAYER_RANGE[0],LAYER_RANGE[1]):
            xemb = xs[layer+1].cpu().detach().numpy()
            vsadd = vs[layer].cpu().detach().numpy()
            Wslayer = Ws[layer].cpu().detach().numpy()
            
            print(xemb.shape)
            print(vsadd.shape)
            print(Wslayer.shape)

            #write xyz and mol file in a temp file
            mol_filename = utils.xyz2mol(props)

            #for each atom in the molecule
            for each_atom in range(number_atoms):
                if each_atom in INDICES_TO_LABEL[each_molecule]:
                    if LABEL == True:
                        #run labeller code to output LDA LABEL, gnuplot color LABEL, gnuplot marker LABEL... NEXT UPDATE!!! 
                        fg_key,ldalabel,gnucolor,gnumarker,hexcolor= labeller.LABEL(mol_filename,number_atoms,each_atom,atomic_numbers[each_atom],each_molecule)
                    else:
                        fg_key = ''
                        ldalabel = 0
                        gnucolor = 0
                        gnumarker = 0
                        hexcolor = ''
                    
                    print(','.join(str(x) for x in Wslayer[0][each_atom][:number_atoms-1,:N_FEATURES]))

                    #write a row of information for the csv file, the row will contaion
                    #molecule_number, the atomic number of atom, ae prediciton, N_FEATURES embedding[layer],
                    #same for the N_FEATURES of vs[layer],
                    row_embs  = ','.join(str(x) for x in xemb[0][each_atom][:N_FEATURES]) +','+str(each_molecule) + ','+str(each_atom) + ','+ str(atomic_numbers[each_atom]) +',' + ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
                    row_vs    = ','.join(str(x) for x in vsadd[0][each_atom][:N_FEATURES]) +','+str(each_molecule) +','+str(each_atom) +  ','+ str(atomic_numbers[each_atom]) +','+ ','.join(str(x) for x in positions[each_atom][:3]) + ','+str(ldalabel) + ','+ str(gnumarker)+','+str(gnucolor)+ ','+ str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom][0].detach().numpy())+','+ fg_key + '\n'
                    flatmatrix_Ws = ','.join(str(x) for x in Wslayer[0][each_atom][:number_atoms-1][:N_FEATURES]) + ',' + str(atomic_numbers[each_atom]) + ',' + ','.join(str(x) for x in positions[each_atom][:3])+','+str(ldalabel)+','+str(gnumarker)+','+str(gnucolor)+','+str(hexcolor)+','+str(layer)+','+str(yi[0][each_atom].detach().numpy()) + ',' + fg_key + '\n'

                    save_xs_datafile.write(row_embs)
                    save_vs_datafile.write(row_vs)
                    save_ws_datafile.write(flatmatrix_Ws)
                else:
                    pass


    save_xs_datafile.close()
    save_vs_datafile.close()
    save_ws_datafile.close()


