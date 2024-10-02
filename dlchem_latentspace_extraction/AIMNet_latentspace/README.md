# DLCheM - What is the chemical logic behing graph neural network (GNN) chemistry


AIMNet_latentspace         - code that extracts fully updated embedding feature vectors from pretrained AIMNet model, which are representative of the latent space of the trained GNN model. 
                             NOTE that labels for the QM9 atom-centered functional groups can be obtained by running the analysis on SchNet first, where  in that code the labelling of atoms is already set with extraction of atom embeddings. Copy the labels after extraction atoms with SchNet and paste them in the resulting .csv here, OR incoporate the label of atoms along with extraction here (by copying the label code and incorporate it with embedding extraction, or just use that code separately on QM9 and copying the resulting label)

                            Inputs required:                     
                                QM9_FILEPATH            - filepath where qm9 db database is being kept
                                QM9_RANGE               - range of QM9 data you want to extract embeddings of
                                AFV_N_FEATURES          - number of features for atomic feature vectors of the pretrained AIMNet model
                                AEF_N_FEATURES          - number of features for atomic environment vectors of the pretrained AIMNet model
                                load_AIMNetMT_ens       - loads trained AIMNet model on wB97x/def2-TZVPP energies, atomic electric moments (charges, dipoles, etc) and volumes
                                load_AIMNetSMD_ens      - loads trained AIMNet model obtained by transfer learning towards SMD-wB97x/def2-TZVPP energies




Requirements: (please check to ensure all dependencies are met! make a conda environment separate for these dependencies if necessary)
    
    python 3.7
    pip install schnetpack==0.3
    conda install pytorch==1.4.0 torchvision==0.5.0 -c pytorch
    pip install numpy==1.21.5
    pip install pandas==1.1.5 
    pip install ipython
    pip install ipykernel
    ipython kernel install --user --name=env_name





