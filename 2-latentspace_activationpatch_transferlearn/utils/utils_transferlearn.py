from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
import math
import joblib

class mlp_train_and_test():
    '''
    mlp_train_and_test ---- class that builds an mlp transfer model from GNN latent space

    3 functions designed to initialize an mlp model, train an mlp model, test an mlp model

    1) init    - set up the scaling fit which uses the training data (must be initialized every time even for testing data to scale on the same fit)
    2) train   - train the model using mlp, and hyperparameters chosen on the training data
    3) test    - test the model on the training dataset 
    '''


    def __init__(self,DATA_FILEPATH,X_COLS_IDXS,Y_COLS_IDXS):
        '''
        Runs the tool to train sklearn's feedforward neural network, with adjustable parameters


        ARGS:
            DATA_FILEPATH           - where the x and y datafile is located, in one file, where the x and y columns will be provided by column indices              
            X_COLS_IDXS             - the column indices of the x data
            Y_COLS_IDXS             - the column indices of the y data
        
        Returns:  
            Standard scaler, fitted to the training data, so that all testing data is scaled by the same scaler fit! 
        '''
        #ARGS
        self.DATA_FILEPATH = DATA_FILEPATH
        self.X_COLS_IDXS = X_COLS_IDXS
        self.Y_COLS_IDXS = Y_COLS_IDXS

        #train_data ---- loading the full dataset as a pandas dataframe   
        train_data = pd.read_csv(self.DATA_FILEPATH)

        #X ---- isolating the part of the dataframe that is for X, using X's column indices specified
        self.X = train_data.iloc[:,X_COLS_IDXS[0]:X_COLS_IDXS[1]].values
        #Y ---- isolating the part of the dataframe that is for Y, using Y's column indices specified
        self.Y = train_data.iloc[:,Y_COLS_IDXS[0]:Y_COLS_IDXS[1]].values

        #scaler_X ---- initializing standard scalar from sklearn for the X data
        self.scaler_X = StandardScaler()
        #scaler_Y ---- initializing standard scalar from sklear for the Y data
        self.scaler_Y = StandardScaler()

        #Check if X has more than one column, or just one column, if one-column, it has to be reshaped for scaling
        if self.X.shape[1] > 1:
            #X_train_scaled ---- normalized X training data using standard scalar on X_train
            self.X_train_scaled = self.scaler_X.fit_transform(self.X)

        #if X has just one column, then it needs to be reshaped for standard scalar because that column would be a numpy row, 
        #and it should be a column with each datapoint as a row
        elif self.X.shape[1] == 1:
            #X_train_scaled ---- normalized X training data using standard scalar on X_train which has ONE col
            self. X_train_scaled = self.scaler_X.fit_transform(self.X.reshape(-1,1))

        #Check if Y has more than one column, or just one column, if one-column, it has to be reshaped
        if self.Y.shape[1] > 1:
            #Y_train_scaled ---- normalized Y training data using standard scalar on Y_train
            self.Y_train_scaled = self.scaler_Y.fit_transform(self.Y)
        elif self.Y.shape[1] == 1:
            self.Y_train_scaled = self.scaler_Y.fit_transform(self.Y.reshape(-1, 1))  


    def train(self,HIDDEN_LAYER_SIZES,ACTIVATION,SAVE_MODEL_FILEPATH):
        '''
        Train mlp on the data 

            HIDDEN_LAYER_SIZES      - tuple to control the number of layers and sizes of the layers (ex. (200,200))
            ACTIVATION              - which activation function to use between the layers
            USE_TRAIN_FOR_TEST      - boolean to choose if to use the training/testing dataset for later testing below,
                                    to compare training vs generalization error 
            SAVE_MODEL_FILEPATH     - filepath here to save the mlp model for later testing

        '''

        #initialize random seed
        seed = 142
        random = np.random.seed(seed)

        #model ---- initializing MLPRegressor model from sklearn with the specified hyperparameters
        model = MLPRegressor(hidden_layer_sizes=HIDDEN_LAYER_SIZES, activation=ACTIVATION, random_state=random,early_stopping=True,validation_fraction=0.1)

        #model ---- running MLPRegressor's fitting attribute to initialized model using the X, Y, scaled training data
        model.fit(self.X_train_scaled, self.Y_train_scaled)

        #predictions_scaled ---- get the normalized predictions by running the fitted model on scaled X_train
        predictions_scaled = model.predict(self.X_train_scaled)

        # Initialize a list to hold MSE for each column
        mse_columns = []

        #Check if Y has more than one column, or just one column, if one-column, predictions have to be reshaped
        if self.Y.shape[1] > 1:
            #predictions ---- inverse scale predictions using standard scaler for Y
            predictions = self.scaler_Y.inverse_transform(predictions_scaled)

            # Calculate MSE for each column
            for i in range(self.Y.shape[1]):
                mse_col = mean_squared_error(self.Y[:, i], predictions[:, i])
                mse_columns.append(mse_col)

            mse_total = mean_squared_error(self.Y, predictions)
            r2 = r2_score(self.Y, predictions)

        elif self.Y.shape[1] == 1:
            #predictions ---- inverse scale predictions using standard scaler for Y, where Y having only one column has to be reshaped
            predictions = self.scaler_Y.inverse_transform(predictions_scaled.reshape(-1,1))

            mse_total = mean_squared_error(self.Y, predictions)
            mse_columns.append(mse_total)  # Since there's only one column, total MSE is the same as the column MSE
            r2 = r2_score(self.Y, predictions)

        #get the validation scores
        mse_model_validation_scores = model.validation_scores_

        #Print results
        print('------------------------------------------------------------------------------')
        print("Model Parameters:")
        print("Hidden layer sizes:", model.hidden_layer_sizes)
        print("Activation function:", model.activation)
        print("Solver:", model.solver)
        print("Learning rate:", model.learning_rate)
        print("Number of iterations:", model.n_iter_)
        print("Coefficients:", model.coefs_)
        print("Number Coefficients:",np.shape(model.coefs_[0]))
        print("Intercepts:", model.intercepts_)
        print("Number Intercepts:",np.shape(model.intercepts_[0]))
        print("Training loss:", model.loss_)
        print("Validation loss scores:", mse_model_validation_scores)

        # Print MSE for each column
        for idx, mse_col in enumerate(mse_columns):
            print(f"RMSE on training data for column {idx + 1}: {math.sqrt(mse_col)}")

        # Print total MSE
        print("Total RMSE on all training data:", math.sqrt(mse_total))
        print("R^2 on run of all training data:",r2)
        print('------------------------------------------------------------------------------')

        #SAVE model
        joblib.dump(model,SAVE_MODEL_FILEPATH)


    def test(self,TEST_DATA_FILEPATH,MODEL_FILEPATH,CENTER_INDICES,CENTERS,NUMBER_CENTERS,COL_TO_ISOLATE,Y_EXISTS,KEY,LABEL_EXISTS,LABEL_COL_INDICES):


        '''
        Testdata preprocessing ---- scaling the testdata for prediction, and inverse scaling predictions to compare with true


        '''

        #test_data ---- loading test data as a dataframe
        test_data = pd.read_csv(TEST_DATA_FILEPATH)


        #if there are no scpeficied center indices, then just get the emebddings of all the atoms in the molecule
        if CENTER_INDICES == None:
            if not COL_TO_ISOLATE:
                #X_test ---- isolating the part of the dataframe that is for X, using X's column indices specified
                X_test = test_data.iloc[:,self.X_COLS_IDXS[0]:self.X_COLS_IDXS[1]].values
                #Y_test ---- isolating the part of the dataframe that is for Y, using Y's column indices specified
                if Y_EXISTS:
                    Y_test = test_data.iloc[:,self.Y_COLS_IDXS[0]:self.Y_COLS_IDXS[1]].values
            else:

                #X_test ---- isolating the part of the dataframe that is for X, using X's column indices specified
                X_test = test_data[test_data[COL_TO_ISOLATE]==KEY].iloc[:,self.X_COLS_IDXS[0]:self.X_COLS_IDXS[1]].values
                
                #Y_test ---- isolating the part of the dataframe that is for Y, using Y's column indices specified
                if Y_EXISTS:
                    Y_test = test_data.iloc[:,self.Y_COLS_IDXS[0]:self.Y_COLS_IDXS[1]].values
        else:
            #indices ---- this is important in case you need to isolate specific indicies for each molecule
            indices = pd.read_csv(CENTER_INDICES,header=None,index_col=None)

            #indices ---- center on whichever indices you want, maybe you do not want all of them
            indices = indices.iloc[:,CENTERS[0]:CENTERS[1]].values


            #X_test ---- initialize the array that ill hold all the x values (embs) for the reaction centers
            X_test = np.zeros((1,128))

            #mol_idx ---- running through each molecule in the number of molecules given given
            for mol_idx in range(NUMBER_CENTERS):
                
                #idx_value ---- getting the molecules's center indices that we want to focus on 
                idx_value = indices[mol_idx]

                #X_idx ---- getting the embedding of the centers of the molecule were interested in
                X_idx = test_data[test_data['mol_idx'] == int(mol_idx)].iloc[int(idx_value),self.X_COLS_IDXS[0]:self.X_COLS_IDXS[1]].values
                
                #X_test ---- stacking the grabbed embeddings to our X embs list to put them all together in one
                X_test = np.vstack((X_test,X_idx))

            #X_test ---- deleting the first row which is just initial zeros
            X_test = np.delete(X_test,0,0) 
        

        #X_test_scaled ---- scaling the X_test, using established scaler for X, so that it is ready for model input
        X_test_scaled = self.scaler_X.transform(X_test)

        #model_loaded --- Load the model from the joblib file it was saved in
        model_loaded = joblib.load(MODEL_FILEPATH)


        #predictions_scaled ---- running the fitted model on the scaled X_test 
        predictions_scaled = model_loaded.predict(X_test_scaled)

        # Check the dimensionality of predictions_scaled
        if len(predictions_scaled.shape) == 1:
            # Reshape 1D array to 2D array with one column
            predictions_scaled = predictions_scaled.reshape(-1, 1)

        #predictions ---- inverse scale predictions using standard scaler for Y
        predictions = self.scaler_Y.inverse_transform(predictions_scaled)

        #Check if we have the Y for this test data, or it is an unknown...
        if Y_EXISTS:
            #mse ---- mean squared error between Y_test and predictions
            mse = mean_squared_error(Y_test, predictions)
            #r2 ---- coefficient of determination between Y_test and predictions
            r2 = r2_score(Y_test, predictions)
            #printing results
            print("RMSE:", math.sqrt(mse))
            print("R^2:",r2)
            print('variance on True:', np.var(Y_test))
            print('Variance on Pred:',np.var(predictions))
        else:
            #printing prediction results
            print('mean of the predictions:',np.mean(predictions))
            print('standard dev of predictions from their mean:', np.std(predictions))

        predictions = pd.DataFrame(predictions)        
        if LABEL_EXISTS and COL_TO_ISOLATE == None:
            label =  test_data.iloc[:,LABEL_COL_INDICES[0]:LABEL_COL_INDICES[1]]
            predictions = pd.concat([predictions,label],axis=1)
        
        return predictions