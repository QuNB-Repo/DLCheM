# DLCheM - What is the chemistry logic behind graph neural network (GNN) chemistry?

How do we analyze "black box" models, specifically graph neural network (GNN) 
models applied to chemistry? Do the patterns they capture after training hold 
any meaningful trends? Can we break down their decision-making? If so, does it 
align with established chemical logic? These are key questions to answer if AI 
is to transform chemical sciences with its speed, efficiency, and reliable 
predictions.

This project explores the latent space of various deep learning neural networks, 
specifically SchNet and AIMNet, to reveal their connection to chemical logic 
and syntax. It examines how this learned logic can be leveraged for transfer 
learning, enabling efficient solutions to new problems in chemistry.

The project consists of three main subprojects:

1) dlchem_latentspace_extraction         
   - This subproject shows how the latent space of a trained GNN model is built 
     on a framework of molecular substructures, also known as functional groups. 
     The framework quantifies molecular similarity between these groups using 
     distance measures. We extract and analyze this latent space using dimensional 
     reduction and distance measuring tools.

2) dlchem_latentspace_reactions          
   - This subproject demonstrates how the latent space follows chemical reaction 
     syntax. Just as natural language models solve word analogies (e.g., "King" - 
     "Man" + "Woman" = "Queen"), a similar approach applies to chemical reactions. 
     For instance, "amides" - "amine" + "alcohol" maps to "carboxylic acid" in the 
     latent space. By constructing reactant and product databases (e.g., from QM9), 
     we analyze reactions using cosine similarity and neighbor tests to uncover 
     chemical reaction analogies within the latent space.

3) dlchem_latentspace_activationpatching_transferlearning   
   - This subproject explores how the latent space, learned from chemical logic, 
     can be used for transfer learning. We apply "activation patching," using 
     activations from a pre-trained model and integrating them into new models 
     to solve problems like pKa, NMR, solubility, and electron density. This 
     approach enables efficient learning with less chemical data and simpler models, 
     which is crucial when reliable data is scarce. It also tests whether deep models 
     can fully capture chemistry, akin to quantum chemistry's wave functions.

Additional projects and endeavors:

4) geneticalgorithms_batteries           
   - Application of genetic algorithms on bipyridine derivatives to predict 
     a bipyridine-like substitute with improved flow battery properties. This 
     study is as an example of successful inverse molecular design.

5) rnn_on_smiles                         
   - 

See the README files within each subproject to learn more!
