# DLCheM - What is the chemical logic behind graph neural network (GNN) chemistry

2) gcnn_latentspace_reactions
   - This subproject demonstrates how the latent space of a trained GNN model 
     follows chemical reaction syntax, similar to how natural language processing 
     (NLP) models capture linguistic patterns. Just as NLP models map word 
     analogies (e.g., "King" - "Man" + "Woman" = "Queen") through vector operations 
     in latent space, we apply the same concept to chemical reactions. For example, 
     the oxidation of alcohol to carbonyl can be represented by a single vector 
     (the oxidation vector, -H2) that applies to all alcohol reactants. This shows 
     that the latent space organizes reactions systematically, allowing us to map 
     reactants to products with consistent vector operations. The algorithm can run 
     any reaction from a database and analyze how the embeddings represent these 
     chemical transformations.

---

**runrxns** 
   - Code that runs reactions across the entire database. It searches for 
     appropriate reactions using the label code (autolabel). Once an appropriate 
     reaction is found, it changes the functional group to the specified product, 
     geometry optimizes it, and returns the dataset of reactants and products 
     for the chosen reaction.

---

**rxns_latentspace** 
   - Jupyter notebook example demonstrating how to analyze reactions in the 
     embedding space. Required inputs:

     - `rxn_centers`: Choose reaction centers from the dictionary of reaction sites 
       (or define your own).
     - `leaving_atoms`: Choose leaving group atoms from the dictionary of leaving groups 
       (or define your own).
     - `attach_group`: Choose a functional group to attach from the dictionary of 
       functional groups to build.
     - `attach_group_bonds`: Select the corresponding bond for the attaching group.
     - `attach_group_numb_atomsbonds`: Specify the number of additional atoms or bonds 
       for the attaching group.
     - `remove_bonds_from_rxnsite`: Optionally remove a double bond from the reaction site.
     - `add_bonds_to_rxnsite`: Optionally add a double bond to the reaction site.
     - `DATASET_FILEPATH`: Path to the labeled database file (via autolabel).
     - `QM9_BOOL`: Boolean indicating whether QM9.db is being used (loads properties 
       with the QM9 method instead of AtomsData for efficiency).
     - `LABELS_FILEPATH`: Path to the file containing atomic center structure labels 
       at sufficient depth (usually a depth of 5 is sufficient).
     - `DEPTH`: Depth of structure label in the label filepath.
     - `MOLECULES_RANGE`: Range of molecules to scan for reactions.
     - `REACT_FILEPATH`: Path to save the reactants.
     - `REMOVELG_FILEPATH`: Path to save reactants without leaving groups (no optimization).
     - `PROD_FILEPATH`: Path to save products with MMFF94 optimization.
     - `AVAILABLE_PROPERTIES`: Placeholder for available properties when making the 
       database. Even if a property is missing, ensure to remember the filler name so 
       it can be accessed.

   The code is well-annotated and explains how to use it at each step. Refer to the 
   annotated Jupyter notebook for further explanations.

---

**embeddings_latentspace** 
   - Code that extracts fully updated embedding feature vectors from SchNet, representing 
     the latent space of the model. It also labels the code with functional group 
     colors/markers and integers for classification.

---

**label** 
   - Code that labels QM9's atoms based on the chemical environment, a.k.a. functional 
     groups around the atom. There are two ways of doing this:
     
     - A) `autolabel` code: Finds each unique atomic center configuration and gives it 
       a unique label with colors and markers.
     - B) `manuallabel` code: Uses a preset label of colors/markers for atomic 
       neighborhood configurations.

---

**Requirements** 
(Please ensure all dependencies are met. Consider creating a separate Conda environment 
for these dependencies):

- Python 3.7
- `conda install -c rdkit rdkit=2020.09.1.0`
- `conda install -c conda-forge openbabel=3.1.1`
- `conda install -c anaconda scikit-learn==1.0.2`
- `conda install -c conda-forge keras=2.2.5`
- `conda install -c anaconda pip=21.2.2`
- `conda install -c anaconda h5py=3.7.0`
- `pip install schnetpack==0.3`
- `pip install numpy==1.21.5`
- `pip install pandas==1.1.5`
- `pip install ipywidgets==8.1.1`
- `pip install captum==0.6.0`
- `pip install ipython`
- `pip install ipykernel`
- `ipython kernel install --user --name=env_name`

---
