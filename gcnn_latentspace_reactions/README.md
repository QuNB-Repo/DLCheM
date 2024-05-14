# DLCheM - What is the chemical logic behing graph neural network (GNN) chemistry

2) gcnn_latentspace_reactions               - This subproject shows how the latent space of a trained GNN model obeys 
                                            chemical reaction syntax, here we run reactions in a database to see 
                                            what would the reactions look like in the latent space of a trained 
                                            GNN model. This reveals that the latent space is organized based on 
                                            reaction syntax, we can drive all reactants of a certain kind to their 
                                            products with a single, constanct vector, for all the reactants. For
                                            example, alcohol to carbonyl oxidation, is represented by a single vector 
                                            for all alcohol reactants (the vector representing oxidation, -H2). 
                                            is true for any reaction chosen. The algorithm in this project is able 
                                            to to run any reactions in a database and analyze the embeddings of

-------------------------------------------------------------------------------------------------------------------

runrxns                                     - code that runs reactions across an entire database. It looks for 
                                            appropriate reactions using the label code (autolabel). Once an 
                                            appropriate reaction is found it changs functional group to specified
                                            product, geometry optimized, and gives back the dataset of reactants and
                                            products for a particular chosen reaction


rxns_latentspace                            - jupyter notebook example of how to analyze reactions in the embedding space
                                            INPUTS required:
                                                rxn_centers          - choosing reaction centers from dictionary of reaction sites (or define own)
                                                leaving_atoms        - choosing leaving groups atoms from dictionary of leaving groups (or define own)
                                                attach_group         - choosing fg to attach from dictionary of fg to build
                                                attach_group_bonds   - choosing the corresponding fg bond to build 
                                                attach_group_numb
                                                _atomsbonds          - the corresponding additional number of atoms/bonds due to the attaching group
                                                remove_bonds_from_
                                                rxnsite              - choosing to remove a double bond from reaction site if required
                                                add_bonds_to_
                                                rxnsite              - choosing to add double bond to reaction site if required
                                                DATASET_FILEPATH     - the db file that has been labelled with autolabel
                                                QM9_BOOL             - a boolean depENDing on whether QM9.db is being used because we
                                                                    load all the properties using QM9 method instead
                                                                    of AtomsData for any other db dataset (more efficient)
                                                LABELS_FILEPATH      - specify the filepath that holds the atomic center structure labels at a sufficient depth
                                                                    for the reaction center (usually a depth of 5 is sufficient for most reaction centers) 
                                                DEPTH                - depth of structure label in label filepath specified
                                                MOLECULES_RANGE      - range of molecules to scan for reactions
                                                REACT_FILEPATH       - filepath that will contain reactants 
                                                REMOVELG_FILEPATH    - filepath that will contain reactants without LG (without optimization)
                                                PROD_FILEPATH        - filepath that will contain products with MMFF94 optimization 
                                                AVAILABLE_PROPERTIES - making the db requires setting an available property even if property is missing
                                                                    this can be anything as filler, long as you remember what the filler database property name is so that it can be accessed

The code is well-annotated and explains how to use it at every step, so feel free to look at the annotated jupyter  notebook for further explanations



--------------------------------------------------------------------------------------------------------------------

embeddings_latentspace     - code that extracts fully updated embedding feature vectors from schnet, 
                             which are representative of the latent space of the model
                             it also labels the code with functional group colors/markers and integers for classification.  

-----------------------------------------------------------------------------------------------------------------------

label                      - code that labels QM9's atoms based on the chemical environment, a.k.a functional groupt 
                             around the atom. There are two ways of doing this: 

                                    A) autolabel code - finds each unique atom center configuraiton and gives it a unique label with colors and markers

                                    B) manuallabel code - preset label of colors/markers for atomic neighborhood configurations 

-----------------------------------------------------------------------------------------------------------------------

Requirements: (please check to ensure all dependencies are met! make a conda environment separate for these dependencies if necessary)
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
    pip install ipywidgets=8.1.1
    pip install captum=0.6.0
    pip install ipython
    pip install ipykernel
    ipython kernel install --user --name=env_name

-----------------------------------------------------------------------------------------------------------------------









			      

