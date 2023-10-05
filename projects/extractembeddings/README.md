**Extract Embeddings README**


To use extractembeddings code you need a datset, we used QM9, which can be downloaded straight from schnetpack using the code we provided (with minor adjustment to allow download=True).

You also need a pre-trained schnet model that has been checkpoint-saved along with its model_state_dict(), it is really important to save the model_state_dict() if you are to train you're own model. 

We have a pre-trained schnet already uploaded and ready for use. However, because it was trained on schnetpack 0.3, it seems to only work on schnetpack 0.3, so that version has to be loaded.
