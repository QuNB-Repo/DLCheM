import pandas as pd
import numpy as np
from tools import utils
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import numpy as np
import math
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error, r2_score
import joblib

np.random.seed(42)

class pertcorrections_neuralnet(nn.Module):
    '''
    simple feedforward neural network which can also reduce to linear regression as well.
    meant to learn only corrections on top of prev perturbational replicates
                        
        Y^{i+1} - Y^{i} =  Y(X_nbrhood_{i})

        we add back the prev perturbational results (Y^{i}) later after the correction is learned

    '''

    def __init__(self,INPUT_DIM,OUTPUT_DIM,NONLINEAR,SECOND_LINEAR,NUMBER_EXTRALAYERS):
        '''
        defines the neural network that serves to add corrections based on nbrhood embedding feature space

            Args:

                INPUT_DIM           - input integer dimension of the neighborhood embedding feature space, this is usually N_feautes*N_depth*5 (for H,C,N,O,F) 
                OUTPUT_DIM          - output integer dimension which is the correction on the previous perturbational embedding
                NONLINEAR           - boolean to decide if to add a non-linearity of tanh to the correction
                SECOND_LINEAR       - boolean to choose another linear layer on top of the linear one above
                NUMBER_EXTRALAYERS  - integer to add extra hidden layers (linear-nonlinear pairs)
            
            returns:
                X^{i+1} while training/testing
        '''


        #super just so we can use this class inside other classes
        super(pertcorrections_neuralnet, self).__init__()

        #setting self.nonlinear to NONLINEAR to be used by forward as well
        self.nonlinear = NONLINEAR

        #same for self.number_extralayers
        self.number_extralayers = NUMBER_EXTRALAYERS

        #setting self.second_linear to SECOND_LINEAR to be used by forward below
        self.second_linear = SECOND_LINEAR

        #self.fc ---- our first layer has to be linear and will serve as outputting the corrections to the prev perturbational embedding
        self.fc = nn.Linear(INPUT_DIM, OUTPUT_DIM)

        #NONLINEAR ---- boolean to decide if too add non-linearity to the corrections
        if self.nonlinear:
            #self.activation ---- adding non-linearity using tanH
            self.activation = nn.Tanh()

        #second_linear ---- if another linear layer is requested after non-linearity, this is always output_dim to output_dim    
        if self.second_linear:
            self.fc2 = nn.Linear(OUTPUT_DIM,OUTPUT_DIM)

        #extra layers of always the same size as the embedding both linear and non-linear
        if self.number_extralayers > 0:
            self.fc3 = nn.Linear(OUTPUT_DIM,OUTPUT_DIM)
    
    def forward(self,x):
        '''
        runs the forward pass for training/testing

            Args:
                x       ---- atom's nbrhood features space input

        NOTE that this only learns the corrections as prev_y has been substracted from Y target
        the prev_y is then added on later as part of the overall algorithm, here we only learn the corrections
        
        '''
        #x ---- linear layer on nbrhood embedding feature space
        x = self.fc(x)

        #x ---- activation on the nbrhood embedding feature space, only if non-linearity is allowed
        if self.nonlinear:
            x = self.activation(x)
        if self.second_linear:
            x = self.fc2(x)
        #x ---- if there are any extra layers (both linear and nonlinear)
        if self.number_extralayers > 0:
            #run extra layers
            for layer_i in range(self.number_extralayers):
                x = self.activation(x)
                x = self.fc3(x)

        #return x for training 
        return x
    
class pertcorrections_neuralnet_regressor(BaseEstimator, RegressorMixin):
    '''
    neural network regression training 

        Args:
            INPUT_DIM           - integer number of input dimensions in the nbrhood embedding features space (Nfeatures*Ndepth*5)
            OUTPUT_DIM          - integer number of output dimensions in the embedding feature space
            NONLINEAR           - boolean to decide if to make the corrections non-linear
            SECOND_LINEAR       - boolean to add another linear on top of the non-linear above
            NUMBER_EXTRALAYERS  - integer to add multiple extra hidden layers (linear-nonlinear)
            EPOCHS              - number of epochs to run the learning algorithm for
            INIT_LR             - learning rate
            BATCHING            - boolean to decide if to batch or not, I found batching useless for this problem
            TRAIN_BATCH_SIZE    - training batch sizes
            VA_BATCH_SIZE       - validation batch sizes
            VALIDATION_SPLIT    - for validation error tracking, ensuring generality acheived during training
            PATIENCE            - for early stopping, epochs to wait if not improvement, then stop
            LR_PATIENCE         - how many epochs to wait before reducing the learning rate
            LR_FACTOR           - factor to reduce learning rate by after waiting on no improvement
            LOG                 - opened log file to write information about training and results
    
    '''
    def __init__(self,INPUT_DIM,OUTPUT_DIM,NONLINEAR,SECOND_LINEAR,NUMBER_EXTRALAYERS,EPOCHS,INIT_LR,BATCHING,TRAIN_BATCH_SIZE,VAL_BATCH_SIZE,VALIDATION_SPLIT,PATIENCE,LR_PATIENCE,LR_FACTOR,LOG):

        #Args
        self.input_dim = INPUT_DIM
        self.output_dim = OUTPUT_DIM
        self.nonlinear = NONLINEAR
        self.second_linear = SECOND_LINEAR
        self.epochs = EPOCHS
        self.lr = INIT_LR
        self.batching = BATCHING
        self.train_batch_size = TRAIN_BATCH_SIZE
        self.val_batch_size = VAL_BATCH_SIZE
        self.validation_split = VALIDATION_SPLIT
        self.patience = PATIENCE
        self.scaler_X = StandardScaler()
        self.scaler_Y = StandardScaler()
        self.model = pertcorrections_neuralnet(INPUT_DIM,OUTPUT_DIM,NONLINEAR,SECOND_LINEAR,NUMBER_EXTRALAYERS)
        self.criterion = nn.MSELoss()
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.lr)
        self.lr_patience = LR_PATIENCE
        self.lr_factor = LR_FACTOR 
        self.scheduler = optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, mode='min', patience=self.lr_patience, factor=self.lr_factor, verbose=True) 
        self.history = {'train_loss': [],
                        'val_loss': []}
        self.log = LOG
        self.number_extralayers = NUMBER_EXTRALAYERS

    def fit(self, X, Y):
        '''
        running the fit between nbrhood embedding feature space and perturbational corrections

            Args:
                X - the neighborhood embedding feature space
                Y - the perturbational correctional targets:  Y_targets - Y_prev_perturbation

        '''

        #X ---- scale X by fitting it to have standard scale (norm of 1)
        X = self.scaler_X.fit_transform(X)
        #Y ---- scale the Y by fitting to have standard scaler 
        # (this is important if you are going to do non-linear!!)
        Y = self.scaler_Y.fit_transform(Y)

        #X ---- X converts to a tensor for pytroch
        X = torch.tensor(X, dtype=torch.float32)
        #Y ---- Y converts to a tensor for pytorch
        Y = torch.tensor(Y, dtype=torch.float32)

        #X_train ---- X data split for the training
        #X_val ---- X data split for the validation
        #Y_train ---- Y data split for the training
        #Y_val ---- Y data split for the validation
        X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size=self.validation_split,random_state=4444)
        #train_dataset ---- training dataset tensor (for batch-loading)
        train_dataset = torch.utils.data.TensorDataset(torch.tensor(X_train,dtype=torch.float32),torch.tensor(Y_train,dtype=torch.float32))
        #val_dataset ---- validation dataset tensor (for batch-loading)
        val_dataset = torch.utils.data.TensorDataset(torch.tensor(X_val, dtype=torch.float32), torch.tensor(Y_val, dtype=torch.float32))
       
        if self.batching == False:
            #train_loader ---- training loader using batches
            train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=len(X_train), shuffle=True)
            #val_loader ---- validation loader using batches 
            val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=len(Y_train), shuffle=True)
        else:
            #train_loader ---- training loader using batches
            train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=self.train_batch_size, shuffle=True)
            #val_loader ---- validation loader using batches 
            val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=self.val_batch_size, shuffle=True)

        #best_val_loss ---- initially a big value
        best_val_loss = float('inf')
        #epochs_no_improve ---- keeping track of the patience
        epochs_no_improve = 0

        #running through the epochs
        for epoch in range(self.epochs):
            #setting model to training mode
            self.model.train()
            #train_loss ---- initially at zero for each epoch
            train_loss = 0.0
            #batch_X ---- extracting randomly shuffled batch of X from train loader
            #batch_y ---- extracting corresponding shuffled batch of Y from the train loader
            for batch_X, batch_Y in train_loader:
                #zero gradients
                self.optimizer.zero_grad()
                #outputs ---- running batches through model
                outputs = self.model(batch_X)
                #loss ---- computing average loss for entire batch
                loss = self.criterion(outputs, batch_Y)
                #backprop gradients through the layers
                loss.backward()
                #updates weights and biases of network using gradients
                self.optimizer.step()
                #train_loss ---- calculating total training loss for entire batch by multiplying average loss by size of the batch
                train_loss += loss.item() * batch_X.size(0)

            #train_loss ---- total loss for added for all items in all batches, now divided by number of datapoints in train_loader dataset
            train_loss /= len(train_loader.dataset)
            #record epoch's train loss
            self.history['train_loss'].append(train_loss)

            #set model to evaluate mode
            self.model.eval()
            #validation loss is initially at zero for the epoch
            val_loss = 0.0
            #without gradient updates
            with torch.no_grad():
                #batch_X ---- extracting randomly shuffled batch of X from val loader
                #batch_Y ----extracting corresponding shuffled batch of X from the val loader
                for batch_X, batch_Y in val_loader:
                    #outputs ---- running batch through model 
                    outputs = self.model(batch_X)
                    #loss ---- average loss on the batch
                    loss = self.criterion(outputs, batch_Y)
                    #val_loss --- calculating total validation loss for the entire batch
                    val_loss += loss.item() * batch_X.size(0)

            #val_loss ---- total loss for added for all items in all batches, now divided by number of datapoints in validation dataset
            val_loss /= len(val_loader.dataset)
            #record epoch's validation loss
            self.history['val_loss'].append(val_loss)

            #update scheduler based on validation loss
            self.scheduler.step(val_loss)

            #print statement every 10 epochs about the training and validation loss
            if (epoch + 1) % 10 == 0 or epoch == self.epochs - 1:
                print(f"Epoch {epoch + 1}/{self.epochs}, Train Loss: {train_loss:.4f}, Validation Loss: {val_loss:.4f}")

            # Check for early stopping
            #if val_loss has improved
            if val_loss < best_val_loss:
                #new record!
                best_val_loss = val_loss
                #start over with the patience count
                epochs_no_improve = 0
            #if val_loss has not improved
            else:
                #add to the patience
                epochs_no_improve += 1

            #if we surpass the patience (number of epochs) STOP training
            if epochs_no_improve >= self.patience:
                print(f"Early stopping at epoch {epoch + 1}")
                #record training/validation results)    
                self.log.write(f'\n          FINAL_TRAINING_RMSE_LOSS:{math.sqrt(train_loss)},     BEST_VAL_RMSE_LOSS:{math.sqrt(best_val_loss)},     NUMBER_EPOCHS: {epoch+1} \n \n')
                self.log.flush()
                break
            elif epoch == self.epochs-1:
                #record training/validation results
                self.log.write(f'\n          FINAL_TRAINING_RMSE_LOSS:{math.sqrt(train_loss)},     BEST_VAL_RMSE_LOSS:{math.sqrt(best_val_loss)},     NUMBER_EPOCHS: {epoch+1} \n \n')               
                self.log.flush()
                
    #running trained model on X for predictions   
    def predict(self, X):
        '''
        Important for subsequent testing of the model 

            Args: 
                X ---- neighborhood embedding feature space
        
            Returns:
                Y_corr ~ Y_target-Y_prev, perturbational corrections 
                
        '''
        #X ---- use fitted standard scaler to scale feature space
        X = self.scaler_X.transform(X)
        #X ---- X converted to tensor
        X = torch.tensor(X, dtype=torch.float32)
        #set model on evaulate mode
        self.model.eval()
        #without updated gradients
        with torch.no_grad():
            #predictions ---- run X through model
            predictions = self.model(X)
            #predictions ---- reverse scale predictions
            predictions = self.scaler_Y.inverse_transform(predictions)
        
        #return predictions
        return predictions
    
    def save_model(self, path):
        #save the model state dictionary so it can be loaded for training/testing later
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'scaler_X': self.scaler_X,
            'scaler_Y': self.scaler_Y,
            'input_dim': self.input_dim,
            'epochs': self.epochs,
            'lr': self.lr,
            'train_batch_size': self.train_batch_size,
            'val_batch_size': self.val_batch_size,
            'validation_split': self.validation_split,
            'patience': self.patience,
            'lr_patience': self.lr_patience,
            'lr_factor': self.lr_factor,
            'history': self.history
        }, path)

    def load_model(self, path):
        #load the model state dictionary if you want to test/train a pretrained model
        checkpoint = torch.load(path)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.scaler_X = checkpoint['scaler_X']
        self.scaler_Y = checkpoint['scaler_Y']
        self.input_dim = checkpoint['input_dim']
        self.epochs = checkpoint['epochs']
        self.lr = checkpoint['lr']
        self.train_batch_size = checkpoint['train_batch_size']
        self.val_batch_size = checkpoint['val_batch_size']
        self.validation_split = checkpoint['validation_split']
        self.patience = checkpoint['patience']
        self.lr_patience = checkpoint['lr_patience']
        self.lr_factor = checkpoint['lr_factor']
        self.history = checkpoint['history']

def constructandfit_nbrlayeremb_vs_GNNemb(LOG_FILE,MOL_RANGE,N_FEATURES,DEPTH,PERTS_RANGE,OG_EMBS_FILEPATH,PERTS_FILEDIR,INPUT_DIM,OUTPUT_DIM,NONLINEAR,SECOND_LINEAR,NUMBER_EXTRALAYERS,EPOCHS,INIT_LR,BATCHING,TRAIN_BATCH_SIZE,VAL_BATCH_SIZE,VALIDATION_SPLIT,PATIENCE,LR_PATIENCE,LR_FACTOR,FIT,PRETRAINED_DIR):
    '''
    constructs nbrhood embedding feature space (with however many depth of layers specified) vs the GNN embedding space for fitting

    the nbrhood embedding feature space must be unique (i.e non-redundant) we order embeddings element-wise and using element-placeholders if an element does not exist
    if an element does not exist as neighbor then it is left as zero for its embedding placeholder, if it exists more than once then its embeddings are added to the same placeholder
    
        Args:
            LOG_FILE             ---- string filepath where to output log info 
            MOL_RANGE            ---- tuple containing start and end range of embedding space to replicate
            N_FEATURES           ---- integer number of features per embedding (dimensionality of embedding space)
            DEPTH                ---- integer specifying how many neighborhood layers deep to use for the neighborhood feature space embedding composition
            PERTS_RANGE          ---- integer specifying number of times to reapply the model, each time applied is based on the previous results,
                                      which leads to incorporating neighborhood information in layers
            OG_EMBS_FILEPATH     ---- string filepath where original embedding vectors located (those to be trained for replication) 
            PERTS_FILEDIR        ---- string filepath where previous perturbational replicates of the training embeddings are located
            INPUT_DIM            ---- number of features per atom's neighborhood feature space
            OUTPUT_DIM           ---- number of features per atom-embedding from GNN model
            NONLINEAR            ---- boolean to decide if to add a non-linear layer on top of the initial linear layer
            SECOND_LINEAR        ---- boolean to decide if to add another linear layer on top of the non-linear
            NUMBER_EXTRALAYERS   ---- integer to add extra hidden layers (linear-nonlinear pairs)
            EPOOCHS              ---- max number of epochs for training
            INIT_LR              ---- initial learning rate
            BATCHING             ---- boolean to decide if to batch the training and validation data (found this useless for this problem, kept this at False)
            TRAIN_BATCH_SIZE     ---- batch size for training
            VAL_BATCH_SIZE       ---- batch size for validation
            VALIDATION_SPLIT     ---- fraction of data to use for validation 
            PATIENCE             ---- how many epochs to wait with no improvement before halting
            LR_PATIENCE          ---- how many epochs to wait with no improvement before reducing the learning rate
            LR_FACTOR            ---- factor by which to reduce the learning rate if no improvement happes over certain patience (epohs) allowed
            FIT                  ---- boolean to decide whether to run a fit or just extract the nbrhood feature space to be used in a different model
            PRETRAINED_DIR       ---- the director of pretrained perturbational models if not fitting but using a pretrained model
    '''

    #embedding_og ---- dataframe of original embedding space
    embedding_og = pd.read_csv(OG_EMBS_FILEPATH,index_col=None)
    
    #Y_nbrhoodembs ---- numpy initial array that will hold each atom's gnn embedding as targets for fitting
    Y_atomembs = embedding_og[(embedding_og['mol_idx'] >= MOL_RANGE[0]) & (embedding_og['mol_idx']  < MOL_RANGE[1])]

    #log ---- Open log file
    log = open(PERTS_FILEDIR+LOG_FILE,mode='a')
    #Record Fitting parameters
    log.write(f'------------------------------------------------Fitting Parameters------------------------------------------------\n')
    log.write(f'          DEPTH: {DEPTH}                       MOL_RANGE: {MOL_RANGE}              N_EMB_FEATURES: {N_FEATURES}  \n\n')
    log.write(f'          DIM_NBRHOOD_EMB: {INPUT_DIM}           DIM_EMB: {OUTPUT_DIM}                         NONLINEAR: {NONLINEAR}  \n\n')
    log.write(f'          MAX_EPOCHS: {EPOCHS}                  INITIAL_LEARNING_RATE: {INIT_LR}           TRAIN_BATCH_SIZE: {TRAIN_BATCH_SIZE}  \n\n')
    log.write(f'          VAL_BATCH_SIZE: {VAL_BATCH_SIZE}            VALIDATION_SPLIT: {VALIDATION_SPLIT}                TOTAL_PATIENCE: {PATIENCE}  \n\n')
    log.write(f'          LR_PATIENCE: {LR_PATIENCE}               LR_FACTOR: {LR_FACTOR}                       NUMBER_ATOMS: {len(Y_atomembs)} \n\n')
    log.write(f'          FIT: {FIT} \n')
    log.write(f'------------------------------------------------Fitting Parameters------------------------------------------------\n')
    log.flush()

    for pert_i in range(PERTS_RANGE[0],PERTS_RANGE[1]):
        
        #record current perturbation
        log.write(f'     PERT: {pert_i} \n')
        log.flush()

        #X_nbrhoodembs ---- numpy initial array that will hold the neighborhood embedding feature space for each perturbation
        X_nbrhoodembs = []

        #Y_pert ---- will hold the final output, 
        #i.e. Y = Y_prev + pert_corr, for each perturbation
        Y_pert = []

        if pert_i !=0:
            prev_Y_filepath = PERTS_FILEDIR+'%s.csv' %(pert_i-1)
            prev_Y = pd.read_csv(prev_Y_filepath,header=None,index_col=False)
            prev_Y = prev_Y[(prev_Y.iloc[:,128] >= MOL_RANGE[0]) & (prev_Y.iloc[:,128]  < MOL_RANGE[1])]


        #molecule_i ---- integer running through all the molecules of the original embedding space, to build the neighborhood composition feature space for each atom in the molecule and return that vs the atom's embedding 
        for molecule_i in range(MOL_RANGE[0],MOL_RANGE[1]):
                
            #simple loading bar
            if molecule_i % 100 == 0:
                print(molecule_i) 

            #embedding_mol_og ---- dataframe isolating each molecules embeddings 
            embedding_mol_og = embedding_og[embedding_og['mol_idx']==molecule_i]

            if pert_i != 0: 
                #prev_y_mol ---- The previous perturbation results for this molecule
                prev_Y_mol = prev_Y[prev_Y.iloc[:,128] == molecule_i].values

            #load some properties for reconstructing xyz of mol, labeled with the og embedding space
            #number_atoms ---- integer number of atoms in the molecule
            number_atoms = len(embedding_mol_og)
            #atomic_numbers ---- dataframe isolate atomic numbers of the atoms in the molecule
            atomic_numbers = embedding_mol_og['atomic_number'].values
            #xyz_positions ---- dataframe isolate xyz coordinates of the atoms in the molecule
            xyz_positions = embedding_mol_og.iloc[:,N_FEATURES+2:N_FEATURES+5].values
            
            #mol_filepath ---- converting xyz and atomic number information above into .mol file where bond connectivity is known so that we can identify direct neighbors
            mol_filepath = utils.convertxyz_to_molfull(atomic_numbers,xyz_positions,TEMP_FILENAME='pertver1')

            #atom_i ---- integer running through total number of atoms in molecule
            for atom_i in range(number_atoms):

                #if pert_i is zero then each atom just gets its average element embedding filled at the right placeholder 
                if pert_i == 0:

                    #atomic_number ---- integer isolating atomic number of atom
                    atomic_number = atomic_numbers[atom_i]
                    
                    #avg_element_emb ---- computing the average element embedding for this atom number as numpy array
                    avg_element_emb = embedding_og[embedding_og['atomic_number'] == atomic_number].iloc[:,0:128].mean().values

                    #X_nbrhoodembs ---- At pert 0, the nbrhood embdding list is just appended with the avg element embeddings corresponding to each atom, that is the zeroth neighborhood layer
                    Y_pert.append(avg_element_emb.tolist()+[int(molecule_i)])


                if pert_i > 0:
                    # organized using element-wise placeholders
                    #here the there are 5 embedding placeholders for each element in QM9 (H,C,N,O,F) in that order, each the size of the embedding dimension 
                    #and multiplied by number of depth of nbrhood layers to include multiple depth of neighborhood composition
                    #if two atoms in a neighborhood layer are the same their embedding is added to the same placeholder,
                    #if an element does not exist as neighbor then it is left as zero for the element's placeholder in that neighborhood layer
                    nbrhood_embed_feature_space = np.zeros((DEPTH*5*N_FEATURES))

                    #nbr_target_lists ---- each layer of neighborhoods has a target list, a list of atom indices that represent the neighbors,
                    #initially, this is just the target atom index, this will help us find its neighbors and replace the target list 
                    #after the second shell (this one below is the first), there will be the possibility of having lists of lists, 
                    #as having two neighbors means each will have their own set of roots
                    nbr_targets_lists = [[atom_i+1]]


                    #for each depth, need to fill the correct embedding placeholders
                    for depth_i in range(DEPTH):

                        #flat_nbrtarget_list ---- as explained above we will get shells of neighbors, because eventually each neighbor has their own family/neighborhood
                        #here we flatten the lists of lists so that we are can just run through all the targets
                        flat_nbrtarget_list = [item for sublist in nbr_targets_lists for item in sublist]

                        #nbr_target_lists ---- cleared so that we now make room for the next layer of neighbors
                        nbr_targets_lists.clear()

                        #run through each target in the neighbor target list to fill their neighborhood embedding feature space at this depth,
                        #then save their neighbors for the next depth
                        for target_i in range(len(flat_nbrtarget_list)):
                             
                            #find nextdoor neighbors using the mol filepath
                            nbrs_idxs, nbrs_ids, nbrs_bond_nmbrs = utils.find_nextdoor_neighbors(mol_filepath,flat_nbrtarget_list[target_i],number_atoms)
                            #nbr_target_lists ---- append all the neighborhood lists for the next perturbation
                            nbr_targets_lists.append(nbrs_idxs) 

                            #number_neighbors ---- number of neighbors for the target
                            number_neighbors = len(nbrs_idxs)

                            #nbrhood_embed_feature_space ---- the atom's neighborhood feature space at specified depth
                            nbrhood_embed_feature_space = utils.fill_nbrhood_featurespace_minimal(nbrhood_embed_feature_space,nbrs_idxs,nbrs_ids,depth_i,N_FEATURES,prev_Y_mol)
                            
                    #X_nbrhoodembs ---- appending the atom's nbrhood feature space to the list of all atom's 
                    #feature spaces as we go through the training embedding vectors to be replicated     
                    X_nbrhoodembs.append(nbrhood_embed_feature_space.tolist())



        '''
        After running through all the molecules to gather the neighborhood feature representation
        Now we fit the feature space to learn the perturbational corrections
        
        '''
        #if pert_i is zero then we are done, we just save the average embeddings
        # as our starting point for perturbational corrections
        if pert_i == 0:
            #Y_pert ---- converting lists of lists into an array
            Y_pert = np.array(Y_pert)

            #replicates_0_df ---- 0th order replicates as a dataframe
            replicates_0_df = pd.DataFrame(Y_pert)
            #save replicate to csv without header nor index
            replicates_0_df.to_csv(PERTS_FILEDIR+'0.csv',index=False,header=False,sep=',')

        #if pert is greater than zero, we then fit perturbational corrections based on neighborhood embedding feature space set up using the prev perturbation
        if pert_i > 0:

            #save model at this perturbation
            pert_filename = '%s' %(pert_i)
            #X_nbrhoodembs_df ---- converting X_nbrhoodembs to dataframe to save to csv
            X_nbrhoodembs_df = pd.DataFrame(X_nbrhoodembs)

            #saving X_nbrhood embs to csv in the perts directory, 
            #this is important incase you want to load the model and run on an input feature space
            X_nbrhoodembs_df.to_csv(PERTS_FILEDIR+pert_filename+'Xnbrhood.csv',index=False,header=False,sep=',')

            #FIT ---- sometimes you want to fit, sometimes you want to use a pretrained model for the updates 
            if FIT == True:

                #Y_corr ---- the perturbational correction to train on
                #need to isolate the original embeddings to match those of X_nbrhoodembs extracted 
                Y_corr = Y_atomembs.iloc[:np.shape(X_nbrhoodembs)[0],0:128].values - prev_Y.iloc[:np.shape(X_nbrhoodembs)[0],0:128].values

    
                #regressor ---- fitting our neural network model
                nn_regressor = pertcorrections_neuralnet_regressor(INPUT_DIM,OUTPUT_DIM,NONLINEAR,SECOND_LINEAR,NUMBER_EXTRALAYERS,EPOCHS,INIT_LR,BATCHING,TRAIN_BATCH_SIZE,VAL_BATCH_SIZE,VALIDATION_SPLIT,PATIENCE,LR_PATIENCE,LR_FACTOR,log)

                #mlp_regressor ---- fitting with mlp
    #            mlp1_regressor = MLPRegressor(hidden_layer_sizes=(OUTPUT_DIM), activation='tanh', random_state=4444)  
    #            mlp2_regressor = MLPRegressor(hidden_layer_sizes=(INPUT_DIM), activation='tanh', random_state=4444)  
    #            mlp3_regressor = MLPRegressor(hidden_layer_sizes=(INPUT_DIM,OUTPUT_DIM), activation='tanh', random_state=4444)  
                                        
                #run fit, return lowest validation error
                nn_regressor.fit(X_nbrhoodembs, Y_corr)
    #            mlp1_regressor.fit(X_nbrhoodembs,Y_corr)
    #            mlp2_regressor.fit(X_nbrhoodembs,Y_corr)
    #            mlp3_regressor.fit(X_nbrhoodembs,Y_corr)


                #save the model in case you want to run it later, using above saved nbrhood feature space
                nn_regressor.save_model(PERTS_FILEDIR+pert_filename+'.pth')
    #            joblib.dump(mlp1_regressor, PERTS_FILEDIR+pert_filename+'mlpver1.pkl')
    #            joblib.dump(mlp2_regressor, PERTS_FILEDIR+pert_filename+'mlpver2.pkl')
    #            joblib.dump(mlp3_regressor, PERTS_FILEDIR+pert_filename+'mlpver3.pkl')

                #corrections ---- get all the corrections at this perturbation after running the all of X through the trained model
                nn_corrections = nn_regressor.predict(X_nbrhoodembs)
    #            mlp1_corrections = mlp1_regressor.predict(X_nbrhoodembs)
    #            mlp2_corrections = mlp2_regressor.predict(X_nbrhoodembs)
    #            mlp3_corrections = mlp3_regressor.predict(X_nbrhoodembs)

                # Evaluate the model
                nn_mse = mean_squared_error(Y_corr,nn_corrections)
    #            mlp1_mse = mean_squared_error(Y_corr, mlp1_corrections)
    #            mlp2_mse = mean_squared_error(Y_corr, mlp2_corrections)
    #            mlp3_mse = mean_squared_error(Y_corr, mlp3_corrections)
                            
                print('nn_rmse: ',math.sqrt(nn_mse))
    #            print('mlp1_rmse: ',math.sqrt(mlp1_mse))
    #            print('mlp2_rmse: ',math.sqrt(mlp2_mse))
    #            print('mlp3_rmse: ',math.sqrt(mlp3_mse))
                log.write(f'nn_rmse:     {math.sqrt(nn_mse)}\n')
    #            log.write(f'mlp1_rmse:    {math.sqrt(mlp1_mse)}\n')
    #            log.write(f'mlp2_rmse:    {math.sqrt(mlp2_mse)}\n')
    #            log.write(f'mlp3_rmse:    {math.sqrt(mlp3_mse)}\n')
                log.write(f'NOTE that discrepancies between final rmse and validation/training rmse \n are probably due to scaling during training and batch-size effects \n Also keep in mind, that with increasing perturbation, we are training on the corrections here, which are getting smaller and smaller and difficult to obtain, \n and therefore the standardized predicted corrections at higher perts will never get that much lower in error compared to standardized correction targets \n O NOT USE IT AS A MEANS TO EVALUATE IF PERTURBATIONS ARE GETTING BETTER AS IT IS SCALED FOR EACH CORRECTION \n\n')


                #replicates ---- obtain replicates by adding corrections to previous perturbational replicate
                replicates = nn_corrections + prev_Y.iloc[:np.shape(X_nbrhoodembs)[0],0:128].values
                #replicates ---- vertical stacking the mol idx of the replicates to be used by further perturbations
                replicates = np.hstack((replicates,Y_atomembs.iloc[:np.shape(X_nbrhoodembs)[0],128:].values))
                #replicates_df ---- converting to pandas dataframe
                replicates_df = pd.DataFrame(replicates)
                #save replicate
                replicates_df.to_csv(PERTS_FILEDIR+pert_filename+'.csv',index=False,header=False,sep=',')
            
            #use a pretrained model to fit the neighborhood feature space to the GNN embeddings
            elif FIT == False:
                
                #Y_corr ---- the perturbational correction to test on
                #need to isolate the original embeddings to match those of X_nbrhoodembs extracted 
                Y_corr = Y_atomembs.iloc[:np.shape(X_nbrhoodembs)[0],0:128].values - prev_Y.iloc[:np.shape(X_nbrhoodembs)[0],0:128].values

                #load_regressor ---- load regressor which will process inputs to outputs, need to load it with all the parameters of the model known
                load_regressor = pertcorrections_neuralnet_regressor(INPUT_DIM,OUTPUT_DIM,NONLINEAR,SECOND_LINEAR,NUMBER_EXTRALAYERS,EPOCHS,INIT_LR,BATCHING,TRAIN_BATCH_SIZE,VAL_BATCH_SIZE,VALIDATION_SPLIT,PATIENCE,LR_PATIENCE,LR_FACTOR,log)
                
                #model_name ---- the name of the model is usually just the perturbation number + .pth
                model_name = '%s.pth' %(pert_i)
                #model_filepath ---- model filepath is where the PRETRAINED model is located plus this current perturbation's name
                model_filepath = PRETRAINED_DIR + model_name
                
                #load_regressor ---- load regressor from pretrained model
                load_regressor.load_model(model_filepath)
                
                #nn_corrections ---- use regressor to get the perturbational updates for the neighborhood
                nn_corrections = load_regressor.predict(X_nbrhoodembs)

                # Evaluate the model
                nn_mse = mean_squared_error(Y_corr,nn_corrections)

                print('nn_rmse: ',math.sqrt(nn_mse))
    #            print('mlp1_rmse: ',math.sqrt(mlp1_mse))
    #            print('mlp2_rmse: ',math.sqrt(mlp2_mse))
    #            print('mlp3_rmse: ',math.sqrt(mlp3_mse))
                log.write(f'nn_rmse:     {math.sqrt(nn_mse)}\n')
    #            log.write(f'mlp1_rmse:    {math.sqrt(mlp1_mse)}\n')
    #            log.write(f'mlp2_rmse:    {math.sqrt(mlp2_mse)}\n')
    #            log.write(f'mlp3_rmse:    {math.sqrt(mlp3_mse)}\n')

                #replicates ---- obtain replicates by adding corrections to previous perturbational replicate
                replicates = nn_corrections + prev_Y.iloc[:np.shape(X_nbrhoodembs)[0],0:128].values
                #replicates ---- vertical stacking the mol idx of the replicates to be used by further perturbations
                replicates = np.hstack((replicates,Y_atomembs.iloc[:np.shape(X_nbrhoodembs)[0],128:].values))
                #replicates_df ---- converting to pandas dataframe
                replicates_df = pd.DataFrame(replicates)
                #save replicate
                replicates_df.to_csv(PERTS_FILEDIR+pert_filename+'.csv',index=False,header=False,sep=',')
            

    #close log file
    log.close()
