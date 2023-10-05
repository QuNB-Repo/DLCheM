**Extract Embeddings README**


To use extractembeddings code you need a datset, we used QM9, which can be downloaded straight from schnetpack using the code we provided (with minor adjustment to allow download, download=True).

You also need a pre-trained schnet model that has been **checkpoint-saved** along with its **model_state_dict()**, it is really important to save the model_state_dict() if you are to train your own model. 

We have a pre-trained schnet already (and model_state_dict() saved), uploaded and ready for use. However, because it was trained on schnetpack 0.3, it seems to only work on schnetpack 0.3, so that version has to be loaded.

Below are the steps to load trained model's state dictionary, extract a layer, such as embedding layers or interaction layers, and perform the same analysis (dimension reduction and lineear discriminant analysis with functional-group-labelled embeddings)

_Step 1) load trained model_

It is really important here to set the exact same hyperparameters used in training the model to load the same model.

You also need to define the dataset filepath (QM9), and load it with schnetpack, if already downloaded, then download = False makes it faster. 

If you saved your model (and model_state_dict()) as we did in our uploaded one, you should have a checkpoint file, using torch to load that checkpoint file, then
load the model state dictionary from the checkpoint

the model state dictionary is printed so you can see the names of all the layers that can be extracted from the trained model 

_Step 2) Define hook function for model layer extraction_

This hook function is a standard way toextract layers from a loaded model. Torch 
has the function "register_forward_hook(hook)." This function registers the defined hook function on a layer (using the layer names from the model state dict), 
so that when the model is applied on an input, the hooked layer can be extracted for that input. 

_Step 3) Extract embeddings from each layer in schnet_

This step sets up a pandas dataframe output that will have the following headings:

#embs0, #embs1, #embs2, #embs3, ...., #embs127, molecule_index, element, x_coord, y_coord, z_coord, fg_label, fg_gnuplotmarker, fg_decimalcolor, layer, fg_key, fg_hexcolor

This part of the code runs the qm9 dataset (but can be replaced with any db dataset using AtomsData from SchNet) through the model, while registering forward hook on initial embedding and interaction layers. The embedding layers can be calculated as emb1 = emb0 + int, where emb0 begins with the initiall "embedding" layer and int are the outputs of each "interaction" layer. 

It saves both the total embedding built (and intermediatte embedding), and the interaction residues that make up the final embedding. 


_Step 5 _

Dimension reduction with sklearn PCA/t-SNE

_Step 6_

linear discriminant analysis with sklearn 






