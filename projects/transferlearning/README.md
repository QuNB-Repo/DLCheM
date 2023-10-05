**Transfer Learning From Embedding to Atomic/Molecular Property**

_data/embsNMR_ - contains the C-embeddings of our NMR dataset (first 200 molecules of QM9 using NMRShiftDB2 model: https://nmrshiftdb.nmr.uni-koeln.de/) and target C-NMR in the 'CembsNPMRDNMR.csv' file

   the embeddings are the first 128 columns, embs0-embs127, 
   target NMR is in column 138. 
   the rest of the columns, 129-137, contain: molecule_index, element (C), x 
   coord, y coord, z coord, fg label, fg gnuplotmarker, fg gnuplotcolor, layer 
   (5), and schnet predicted atomization energy 

_data/embspKa_ - contains the O-embeddings of our pKa dataset (curated from IUPAC/pKa database: https://github.com/IUPAC/Dissociation-Constants) and target pKas in the '601embs.csv' file

   the embeddings are the first 128 columns, embs0-embs127
   target pKa is in column 128. 
   the rest of the columns, 129-133 contain the ldalabel, gnuplotmarker, 
   gnuplotcolor of the oxygen atom

_data/embslogS _ - contains the embeddings of 800 molecules along with their logS solubility curated from the NP-MRD database found here: https://np-mrd.org/. 

   Embeddings file
   the embeddings are the first 128 columns, embs0-embs127, 
   the rest of the columns, 129-137, contain: molecule_index, element (C), x 
   coord, y coord, z coord, fg label, fg gnuplotmarker, fg gnuplotcolor, layer 
   (5), and schnet predicted atomization energy 


1) lrregression - Linear regression 

Needs an input embeddings with target property (we put them in the same csv file, found in data/)

   
3) mlpregression - sklearn's multi-layer perceptron (feedforward neural net) 

Needs an input embeddings with target property (we put them in the same csv file, found in data/)

   
4) sizeextnn - size-extensive neural net a neural network design that trains on the sum of atomisitic contributions from all embeddings of the molecule
