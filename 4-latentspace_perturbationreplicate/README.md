# Latent Space Perturbation Replication - Replicating the Latent Space of GNN Models

How is a graph latent space constructed? Can we reconstruct it?

This project demonstrates how to replicate the latent space of graph 
neural network (GNN) models by iteratively fitting perturbational 
replicates. The procedure incrementally updates the node embeddings 
based on neighborhood information and refines the representation 
through multiple iterations.

## Procedure:

1) **Initialize** all nodes (per molecule) in the dataset with a guess 
   (e.g., the mean nodal activation of the dataset).
2) **First-order fit** of the latent space based on neighborhood nodal 
   vectors (using a cutoff distance).
3) **Reapply the fit** of the nodal latent vector based on the results 
   of the previous fit. This step is critical, allowing fine-tuning of 
   the latent space replicates by incorporating deeper atom-neighbor 
   and neighbor-neighbor correlations.
4) **Repeat** step 3 until an accurate latent space replicate is 
   achieved.

---

## Project Structure

### 1) **Neighborhood Feature Construction and Fitting Phase**

This phase combines the construction of the neighborhood feature 
embedding space and fitting the perturbational replicate embeddings 
in one method. The reason for combining these tasks is that the 
neighborhood embeddings need to be updated for each perturbation 
step to ensure consistency. 

**Inputs required for neighborhood feature construction:**

- `LOG_FILE`: Filepath to output log information.
- `MOL_RANGE`: Range of molecules (start, end) for which to replicate 
  the embedding space.
- `N_FEATURES`: Number of features per embedding (dimensionality of 
  the embedding space).
- `DEPTH`: Number of neighborhood layers deep to use for the neighborhood 
  feature space embedding composition.
- `PERTS_RANGE`: Number of perturbational iterations to apply, each 
  iteration reuses the result from the previous.
- `OG_EMBS_FILEPATH`: Path to the original embedding vectors for 
  replication.
- `PERTS_FILEDIR`: Directory where perturbational replicate embeddings 
  will be saved.

**Inputs required for fitting:**

- `INPUT_DIM`: Dimensions per atom's neighborhood features (e.g., 
  `NUMBER_ELEMENTS * DEPTH * DIM_EMBEDDING`).
- `OUTPUT_DIM`: Output dimensions per atom's embedding from the 
  pretrained GNN model.
- `NONLINEAR`: Boolean to indicate if a nonlinear layer (tanh) should 
  be added.
- `SECOND_LINEAR`: Boolean to decide if an extra linear layer should 
  be added.
- `NUMBER_EXTRALAYERS`: Number of extra layers to test for hyperparameters.
- `EPOCHS`: Maximum number of epochs for training.
- `INIT_LR`: Initial learning rate for backpropagation.
- `BATCHING`: Boolean to decide whether to batch training and validation 
  (though this was found unnecessary in this study).
- `TRAIN_BATCH_SIZE`: Training batch size if batching is True.
- `VAL_BATCH_SIZE`: Validation batch size if batching is True.
- `VALIDATION_SPLIT`: Fraction of total data to use for validation.
- `PATIENCE`: Number of epochs to wait without improvement before 
  halting training.
- `LR_PATIENCE`: Number of epochs to wait before reducing the learning 
  rate.
- `LR_FACTOR`: Factor by which to reduce the learning rate.
- `FIT`: Boolean to determine if fitting should occur, or if only 
  neighborhood feature representation should be built using a pretrained model.
- `PRETRAINED_DIR`: Directory of pretrained perturbational models if 
  using a pretrained model instead of fitting.

---

### 2) **Fitting/Training Phase**

In this phase, the model is trained to fit the perturbational replicate 
embeddings based on neighborhood features. This process is repeated 
until all perturbations have been applied, progressively refining the 
neighborhood embedding composition and fitting to the GNN embeddings.

---

## Requirements

To set up the environment, the following dependencies are required:

```bash
python 3.7
pandas==1.1.5
numpy==1.21.5
torch==1.9.0
scikit-learn==1.0.2
