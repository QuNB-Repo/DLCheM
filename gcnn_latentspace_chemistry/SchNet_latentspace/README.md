# DLCheM - What is the chemical logic behing graph neural network (GNN) chemistry


SchNet_latentspace         - code that extracts fully updated embedding feature vectors from schnet, which are representative of the latent space of the model
                             it also labels the code with functional group colors/markers and integers for classification.  

                            Inputs required:                     
                                MODEL_FILEPATH              - where the pretrained schnet model is saved
                                DB_FILEPATH                 - where the molecules db file for embedding extraction is saved
                                SAVE_FILEPATH               - where to save the extracted and labelled embeddings
                                QM9_TRUE                    - boolean to load qm9.db faster (otherwise you will have to list all available properties)
                                AVAILABLE_PROPERTIES        - list of the available properties in the db database ['property1','property2']
                                START                       - index of start molecule of embedding extraction from the db
                                END                         - index of the end molecule of embedding extraction from the db 
                                N_FEATURES                  - number of features in embedding (depends on model trained)
                                LAYERS                      - layers (range) to extract (depends on number of layers in model, can be less, not more)
                                ELEMENTS                    - list of elements to extract embeddings for (doesn't have to be all of the existing elements in db)
                                LABEL                       - boolean to decide if to label the embeddings with a FG label colors and markers
                                RESTRICT_LABEL              - boolean to decide if to only extract embeddings of those with a certain FG label (if you know the label prehand)
                                ALLOWED_LABELS              - list of allowed FG labels (integers) to extract embeddings for
                                SCRATCH_FILE                - a scratch file to write xyz and mol onto temporarily for labelling FG's
                                ADD_HEADER                  - boolean to add a header to the output embedding csv file



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



