# Latent Space Extraction of SchNet Model

---

## embeddings_latentspace
This code extracts the updated embedding feature vectors from a pretrained 
SchNet model, representing the latent space of the model. It also labels 
the embeddings with functional group (FG) colors, markers, and integers 
for classification.

---

## extractembs.ipynb
This Jupyter notebook example demonstrates the usage of the 
`embeddings_latentspace` code to extract the latent space embeddings.

**Inputs required:**

- `MODEL_FILEPATH`: Path to the pretrained SchNet model.
- `DB_FILEPATH`: Path to the molecule database file for embedding extraction.
- `SAVE_FILEPATH`: Path to save the extracted and labeled embeddings.
- `QM9_TRUE`: Boolean flag for loading the QM9 database faster. If False, 
  you must specify all available properties.
- `AVAILABLE_PROPERTIES`: List of available properties in the database 
  (e.g., `['property1', 'property2']`).
- `START`: Index of the starting molecule for embedding extraction.
- `END`: Index of the ending molecule for embedding extraction.
- `N_FEATURES`: Number of features in the embedding (depends on the 
  model used for training).
- `LAYERS`: Range of layers to extract (cannot exceed the number of 
  layers in the model).
- `ELEMENTS`: List of elements to extract embeddings for (not necessary 
  to include all elements in the database).
- `LABEL`: Boolean flag to decide whether to label the embeddings with 
  FG labels (colors and markers).
- `RESTRICT_LABEL`: Boolean flag to restrict extraction to embeddings 
  with a certain FG label (if known beforehand).
- `ALLOWED_LABELS`: List of allowed FG labels (integers) to extract 
  embeddings for.
- `SCRATCH_FILE`: A scratch file used temporarily to write XYZ and MOL 
  data for labeling FG's.
- `ADD_HEADER`: Boolean flag to add a header to the output embedding CSV.

---

## label
This code labels QM9 atoms based on their chemical environment (functional 
groups) around the atom. Two methods are available for labeling:

1) `autolabel`: Automatically identifies each unique atomic center 
   configuration and assigns a unique label with colors and markers.
2) `manuallabel`: Uses preset labels with colors and markers for 
   atomic neighborhood configurations.

---

## Requirements

To set up the environment, ensure all dependencies are met. It is 
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
