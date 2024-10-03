# Latent Space Reactions in Deep Learning Models
---

## gcnn_latentspace_reactions
This subproject demonstrates how the latent space of a trained GNN model 
follows chemical reaction syntax, much like how natural language processing 
(NLP) models capture linguistic patterns. Just as NLP models map word 
analogies (e.g., "King" - "Man" + "Woman" = "Queen") through vector operations 
in latent space, we apply this concept to chemical reactions. For example, 
the oxidation of alcohol to carbonyl is represented by a single vector 
(oxidation vector, -H2) that applies to all alcohol reactants. This shows 
how the latent space organizes reactions systematically, allowing reactants 
to be mapped to products through consistent vector operations. The algorithm 
can run any reaction from a database and analyze how the embeddings represent 
these chemical transformations.

---

## runrxns
This module runs reactions across the entire database. It searches for 
appropriate reactions using a label code (autolabel). Once a reaction is 
found, it alters the functional group to the specified product, performs 
geometry optimization, and returns the dataset of reactants and products 
for the selected reaction.

---

## rxns_latentspace
This Jupyter notebook provides an example of how to analyze reactions 
within the embedding space. The following inputs are required:

- `rxn_centers`: Choose reaction centers from a reaction site dictionary 
  (or define your own).
- `leaving_atoms`: Choose leaving group atoms from a dictionary of leaving 
  groups (or define your own).
- `attach_group`: Choose a functional group to attach from the dictionary 
  of functional groups to build.
- `attach_group_bonds`: Select the corresponding bond for the attaching group.
- `attach_group_numb_atomsbonds`: Specify the additional number of atoms or 
  bonds for the attaching group.
- `remove_bonds_from_rxnsite`: Optionally remove a double bond from the 
  reaction site.
- `add_bonds_to_rxnsite`: Optionally add a double bond to the reaction site.
- `DATASET_FILEPATH`: Path to the labeled database file (via autolabel).
- `QM9_BOOL`: Boolean flag indicating if QM9.db is used (loads properties 
  via QM9 method for efficiency).
- `LABELS_FILEPATH`: Path to the atomic center structure labels file.
- `DEPTH`: Depth of structure label in the label filepath.
- `MOLECULES_RANGE`: Range of molecules to scan for reactions.
- `REACT_FILEPATH`: Path to save the reactants.
- `REMOVELG_FILEPATH`: Path to save reactants without leaving groups 
  (without optimization).
- `PROD_FILEPATH`: Path to save products with MMFF94 optimization.
- `AVAILABLE_PROPERTIES`: Placeholder for available properties during 
  database creation. Even if a property is missing, remember its filler name 
  for future access.

The code is well-annotated, explaining each step. Refer to the annotated 
Jupyter notebook for more details.

---

## embeddings_latentspace
This code extracts fully updated embedding feature vectors from a 
pretrained SchNet model, representing the latent space of the model. It 
also labels the embeddings with functional group colors/markers and 
integers for classification.

---

## label
This code labels QM9 atoms based on their chemical environment, or 
functional groups around the atom. Two approaches are available for 
labeling:

- A) `autolabel`: Automatically assigns a unique label, with colors and 
  markers, to each unique atomic center configuration.
- B) `manuallabel`: Uses preset labels with colors and markers for atomic 
  neighborhood configurations.

---

## Requirements

To set up the environment, ensure all dependencies are installed. It is 
recommended to create a separate Conda environment for this project.

```bash
python 3.7
conda install -c rdkit rdkit=2020.09.1.0
conda install -c conda-forge openbabel=3.1.1
conda install -c anaconda scikit-learn==1.0.2
conda install -c conda-forge keras=2.2.5
conda install -c anaconda pip=21.2.2
conda install -c anaconda h5py=3.7.0
pip install schnetpack==0.3
pip install numpy==1.21.5
pip install pandas==1.1.5
pip install ipywidgets==8.1.1
pip install captum==0.6.0
pip install ipython
pip install ipykernel
ipython kernel install --user --name=env_name
