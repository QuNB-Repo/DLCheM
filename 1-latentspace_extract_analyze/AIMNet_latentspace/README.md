# AIMNet Latent Space Extraction

This code extracts the embedding feature vectors from a pretrained 
AIMNet model, representing the latent space of the trained GNN model. 

**Important**: Labels for the QM9 atom-centered functional groups can 
be obtained by first running the extraction process on SchNet, as that 
code includes labeling for atoms during the embedding extraction. After 
extracting atom embeddings with SchNet, you can copy the resulting labels 
into the output CSV from this AIMNet extraction. Alternatively, you can 
integrate the labeling code into this script, or run the label extraction 
separately and copy the labels over to this file.

---

## Inputs required:

- `QM9_FILEPATH`: Path to the QM9 database file.
- `QM9_RANGE`: Range of QM9 data (start, end) to extract embeddings from.
- `AFV_N_FEATURES`: Number of features for atomic feature vectors from 
  the pretrained AIMNet model.
- `AEF_N_FEATURES`: Number of features for atomic environment vectors 
  from the pretrained AIMNet model.
- `load_AIMNetMT_ens`: Loads a pretrained AIMNet model trained on 
  wB97x/def2-TZVPP energies, atomic electric moments (charges, dipoles, 
  etc.), and volumes.
- `load_AIMNetSMD_ens`: Loads a pretrained AIMNet model obtained via 
  transfer learning towards SMD-wB97x/def2-TZVPP energies.

---

## Requirements

To set up the environment, please ensure all dependencies are met. 
It is recommended to create a separate Conda environment for this project.

```bash
python 3.7
pip install schnetpack==0.3
conda install pytorch==1.4.0 torchvision==0.5.0 -c pytorch
pip install numpy==1.21.5
pip install pandas==1.1.5
pip install ipython
pip install ipykernel
ipython kernel install --user --name=env_name
