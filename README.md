# DLCheM - What is the chemistry logic behind deep learning models in chemistry?

How do we analyze "black box" models, specifically deep learning models applied 
to chemistry? Do the patterns they capture after training hold any meaningful 
trends? Can we break down their decision-making? If so, does it align with 
established chemical logic? These are key questions to answer if AI is to 
transform chemical sciences with its speed, efficiency, and reliable 
predictions.

This project explores the latent space of various deep learning neural networks, 
specifically SchNet and AIMNet, to reveal their connection to chemical logic 
and syntax. It examines how this learned logic can be leveraged for transfer 
learning, enabling efficient solutions to new problems in chemistry.

The project consists of three main subprojects:

1) latentspace_extract_analyze         
   - This subproject shows how the latent space of a trained GNN model is built 
     on a framework of molecular substructures, also known as functional groups. 
     The framework quantifies molecular similarity between these groups using 
     distance measures. We extract and analyze this latent space using dimensional 
     reduction and distance measuring tools.

     read README file inside project directory for more details

     publication: A.M. El-Samman, I.A. Husain, M. Huynh, S. De Castro, B. Morton, 
                  and S. De Baerdemacker. Global geometry of chemical graph neural 
                  network representations in terms of chemical moieties. Digital 
                  Discovery, 3(3):544–557, 2024

2) latentspace_reaction_transformer        
   - This subproject demonstrates how the latent space follows chemical reaction 
     syntax. Just as natural language models solve word analogies (e.g., "King" - 
     "Man" + "Woman" = "Queen"), a similar approach applies to chemical reactions. 
     For instance, "amides" - "amine" + "alcohol" maps to "carboxylic acid" in the 
     latent space. By constructing reactant and product databases (e.g., from QM9), 
     we analyze reactions using cosine similarity and neighbor tests to uncover 
     chemical reaction analogies within the latent space.

     read README file inside project directory for more details

     publication:  A. M. El-Samman and S. De Baerdemacker, “amide - amine + alcohol
                   = carboxylic acid. Chemical reactions as linear algebraic analogies 
                   in graph neural networks,” ChemRxiv. 2024; doi:10.26434/chemrxiv-2024-fmck4

3) latentspace_activationpatch_transferlearn  
   - This subproject explores how the latent space, learned from chemical logic, 
     can be used for transfer learning. We apply "activation patching," using 
     activations from a pre-trained model and integrating them into new models 
     to solve problems like pKa, NMR, solubility, and electron density. This 
     approach enables efficient learning with less chemical data and simpler models, 
     which is crucial when reliable data is scarce. It also tests whether deep models 
     can fully capture chemistry, akin to quantum chemistry's wave functions.

     read README file inside project directory for more details

     publication: A. M. El-Samman, S. De Castro, B. Morton and S. De Baerdemacker. 
                  Transfer learning graph representations of molecules for pKa, 
                  13C-NMR, and solubility. Canadian Journal of Chemistry, 2023, 
                  102, 4.

4) latentspace_perturbation_replicate
   - This subproject is a method to replicate the latent space of graph neural networks
     using an input feature space based on local chemical neighborhood composition. 
     if iteratively repeated based on the results of the previous replicate you get
     perturbational corrections to the latent space replicate. The reason this happens
     is because each neighbor is being updated with its neighbors such that with
     greater iterations neighbor-neighbor and atom-neighbor correlations are introduced
     to replicate the latent space.

     full details on this in the upcoming manuscript

     publication: manuscript in progress

     read README file inside project directory for more details


See the README files within each subproject to learn more!

5) Other projects: folder contains other projects and endeavors! Which are:

  geneticalgorithms_batteries           
   - Application of genetic algorithms on bipyridine derivatives to predict 
     a bipyridine-like substitute with improved flow battery properties. This 
     study is as an example of successful inverse molecular design.

  rnn_onsmiles                         
   - a reimplementation of smiles2vec to analyze the latent space of an RNN-type model


  Again, see the README files within each subproject to learn more!
