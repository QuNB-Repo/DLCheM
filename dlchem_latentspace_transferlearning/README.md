# 1) Activation Patching Transfer Learning (see below for Activation Patching SizeExtNN Transfer Learning)

This project implements a transfer learning framework using activation 
patching to solve various chemical property prediction tasks. The model 
is initially trained on a molecular property, and the latent space is 
extracted using the `extractembs.ipynb` project. Once the latent space 
is obtained, this code is used to train a model where the X data 
represents the latent space of the training data, and the Y data 
represents the property to predict.

## Project Structure

### 1) **Training Phase**

The training phase initializes the preprocessing of the latent space 
data, scales the input features, and trains a multi-layer perceptron (MLP) 
model. The trained model is saved for future use in the testing phase.

**Inputs required for training:**

- `DATA_FILEPATH`: Path to the CSV file containing the training latent 
  space data, with both input (X) and target (Y) columns.
- `X_COLS_IDXS`: List of column indices representing the features (X), 
  which are the latent space vectors extracted from the training data.
- `Y_COLS_IDXS`: List of column indices representing the target variables 
  (Y) to predict (e.g., molecular properties like pKa or electron density).
- `SAVE_PATH`: Path to save the trained model for later use.

The code will save a `.joblib` model file after training, which is used 
in the next testing phase.

---

### 2) **Testing Phase**

The testing phase uses the trained model to make predictions on a new 
dataset. To ensure consistency, the original training dataset must be 
loaded to retrieve the scaler function applied during the training phase.

**Important**: Even when testing on a new dataset, the original training 
dataset must be loaded to apply the same scaling function. You will need 
to load the `DATA_FILEPATH` from the training phase along with the 
same `X_COLS_IDXS` and `Y_COLS_IDXS` used during training to correctly 
retrieve the scaler function.

**Inputs required for testing:**

- `TEST_DATA_FILEPATH`: Path to the CSV file containing the testing latent 
  space data.
- `MODEL_FILEPATH`: Path to the saved model from the training phase.
- `Y_EXISTS`: Boolean flag indicating if the target values exist in the 
  testing data (for evaluation purposes).
- `LABEL_EXISTS`: Boolean flag indicating if label columns exist in the 
  testing data.
- `COL_TO_ISOLATE`: Column name or identifier to isolate for analysis 
  during prediction.
- `ROW_KEY`: Row index to focus on within the dataset.
- `LABEL_COL_IDXS`: Column indices for labels if they exist.
- `CENTER_INDICES_FILEPATH`: Optional file path for center indices if 
  required.
- `CENTERS`: Tuple specifying the center points.
- `NUMBER_CENTERS`: Number of centers for the prediction task.

After testing, the predictions are saved to a specified output path.

---

## Latent Space Extraction

Before using this code, the latent space of the training data must be 
obtained using the `extractembs.ipynb` project. This step is crucial, as 
the X data in the training and testing phases represents the extracted 
latent space vectors. These vectors capture important chemical features 
from the training molecules, and the Y data corresponds to the specific 
property being predicted (e.g., pKa, NMR, solubility).

---

## Requirements

To set up the environment, you will need the following dependencies:

```bash
python 3.7
torch==1.9.0
conda install -c anaconda scikit-learn==1.0.2
pip install numpy==1.21.5
pip install pandas==1.1.5


# Size-Extensive Neural Network Activation Patching Transfer Learning

This project implements a size-extensive neural network architecture 
for activation patching transfer learning. The network leverages the 
latent space embeddings of a pretrained model to fine-tune and predict 
chemical properties. This approach sums the atom-wise contributions 
to predict molecular properties, ensuring size-extensivity in the 
network's architecture.

## Project Structure

### 1) **Neural Network Architecture**

The neural network is built with a size-extensive architecture, 
where each atom's latent space embedding is passed through a 
two-layer fully connected neural network. The atom-wise results 
are summed to compute the final predicted molecular property.

**Inputs required for the network:**

- `INPUT_SIZE`: The number of dimensions in the input feature 
  space (latent space embeddings).
- `HIDDEN_SIZE`: The number of parameters in the hidden layer.
- `OUTPUT_SIZE`: The number of output dimensions (typically 1 
  for a scalar property like molecular energy).
  
The network uses ReLU activation in the hidden layer and 
summed outputs for the final prediction.

---

### 2) **Training Phase**

The training phase involves loading the latent space embeddings 
and the molecular property data, normalizing the embeddings, 
and training the network using a mean squared error loss function 
and the Adam optimizer.

**Inputs required for training:**

- `EMBS_PATH`: Path to the CSV file containing the latent 
  space embeddings.
- `N_EMBS_FEATURES`: Number of features in the latent space 
  embeddings (typically 128).
- `MOL_PROPERTY_PATH`: Path to the CSV file containing the 
  molecular properties to be predicted.
- `INPUT_SIZE`: Number of dimensions in the input latent space 
  (128 in this case).
- `HIDDEN_SIZE`: Number of parameters in the hidden layer (200 
  in this case).
- `OUTPUT_SIZE`: Number of output dimensions (1 for scalar property).
- `LEARNING_RATE`: The learning rate for the Adam optimizer 
  (0.001 in this case).
- `NUM_EPOCHS`: Number of epochs for training (10,000 in this case).
- `NUM_TRAIN_SAMPLES`: Number of training samples (450 molecules).
- `NUM_VAL_SAMPLES`: Number of validation samples (100 molecules).

The training loop updates the model’s weights by minimizing 
the loss between predicted and true molecular properties, using 
the Adam optimizer.

---

### 3) **Testing Phase**

After training, the model is evaluated on the validation dataset 
to calculate the validation loss. The latent space embeddings of 
the validation molecules are fed through the model, and the 
predictions are compared against the true molecular properties.

**Important**: The embeddings used for training must be 
normalized, and the same normalized embeddings should be used 
in the testing phase.

**Inputs required for testing:**

- `TEST_EMBS_PATH`: Path to the CSV file containing the testing 
  latent space embeddings.
- `MOL_PROPERTY_PATH`: Path to the CSV file containing the true 
  molecular properties for validation.
- `MODEL_FILEPATH`: Path to the saved model after training.

During validation, the average loss is printed for each epoch 
to monitor performance on unseen data.

---

## Requirements

To set up the environment, the following dependencies are required:

```bash
python 3.7
torch==1.9.0
numpy==1.21.5
pandas==1.1.5


## 2) Size-Extensive Neural Network Activation Patching Transfer Learning

This project implements a size-extensive neural network architecture 
for activation patching transfer learning. The network leverages the 
latent space embeddings of a pretrained model to fine-tune and predict 
chemical properties. This approach sums the atom-wise contributions 
to predict molecular properties, ensuring size-extensivity in the 
network's architecture.

## Project Structure

### 1) **Neural Network Architecture**

The neural network is built with a size-extensive architecture, 
where each atom's latent space embedding is passed through a 
two-layer fully connected neural network. The atom-wise results 
are summed to compute the final predicted molecular property.

**Inputs required for the network:**

- `INPUT_SIZE`: The number of dimensions in the input feature 
  space (latent space embeddings).
- `HIDDEN_SIZE`: The number of parameters in the hidden layer.
- `OUTPUT_SIZE`: The number of output dimensions (typically 1 
  for a scalar property like molecular energy).
  
The network uses ReLU activation in the hidden layer and 
summed outputs for the final prediction.

---

### 2) **Training Phase**

The training phase involves loading the latent space embeddings 
and the molecular property data, normalizing the embeddings, 
and training the network using a mean squared error loss function 
and the Adam optimizer.

**Inputs required for training:**

- `EMBS_PATH`: Path to the CSV file containing the latent 
  space embeddings.
- `N_EMBS_FEATURES`: Number of features in the latent space 
  embeddings (typically 128).
- `MOL_PROPERTY_PATH`: Path to the CSV file containing the 
  molecular properties to be predicted.
- `INPUT_SIZE`: Number of dimensions in the input latent space 
  (128 in this case).
- `HIDDEN_SIZE`: Number of parameters in the hidden layer (200 
  in this case).
- `OUTPUT_SIZE`: Number of output dimensions (1 for scalar property).
- `LEARNING_RATE`: The learning rate for the Adam optimizer 
  (0.001 in this case).
- `NUM_EPOCHS`: Number of epochs for training (10,000 in this case).
- `NUM_TRAIN_SAMPLES`: Number of training samples (450 molecules).
- `NUM_VAL_SAMPLES`: Number of validation samples (100 molecules).

The training loop updates the model’s weights by minimizing 
the loss between predicted and true molecular properties, using 
the Adam optimizer.

---

### 3) **Testing Phase**

After training, the model is evaluated on the validation dataset 
to calculate the validation loss. The latent space embeddings of 
the validation molecules are fed through the model, and the 
predictions are compared against the true molecular properties.

**Important**: The embeddings used for training must be 
normalized, and the same normalized embeddings should be used 
in the testing phase.

**Inputs required for testing:**

- `TEST_EMBS_PATH`: Path to the CSV file containing the testing 
  latent space embeddings.
- `MOL_PROPERTY_PATH`: Path to the CSV file containing the true 
  molecular properties for validation.
- `MODEL_FILEPATH`: Path to the saved model after training.

During validation, the average loss is printed for each epoch 
to monitor performance on unseen data.

---

## Requirements

To set up the environment, the following dependencies are required:

```bash
python 3.7
torch==1.9.0
numpy==1.21.5
pandas==1.1.5

