# DLCheM - What is the chemical logic behing graph neural network (GNN) chemistry

How do we analyze "black box" models, specifically grah neural networks (GNN) models on chemistry? Does the information they store after training hold any meaningful trends? Is it possible to break down the decision making of such models, and if so, is it using proper chemical logic as part of its decision-making framework? These are the questions that must be tackled if AI is to revolutionalize the chemical sciences with its speed, efficiency, and trustworthy accurate predictions.

This project analyzes the latent space of graph convolutional neural networks, in particular SchNet and AIMNet, 
to reveal their connection to chemical logic and chemical syntax, and how this logic learned by the algorithm can be used to transfer solve new problems in chemistry very efficiently! 

The project is divided in three subprojects: 

    1) gcnn_latentspace_chemistry               - This subproject shows how the latent space of a trained model
                                                  is related to a chemical environment decisions-making framework based on molecular environments, a.k.a functional groups. It is so precise, it quantifies the difference between local molecular shapes. This allows the algorithm to identify the local structures of chemistry, thus basing the decision on this fundamental concept. 
    2) gcnn_latentspace_reactions               - This subproject shows how the latent space of a trained model obeys 
                                                  chemical reaction syntax, here we run reactions in a database to see what would the reactions look like in the latent space of a trained GNN model. This reveals that the latent space is organized based on reaction syntax, we can drive all reactants of a certain kind to their products with a single, constanct vector, for all the reactants. This is true for any reaction chosen. The algorithm in this project is able to to run any reactions in a database and analyze the embeddings of
    3) gcnn_latentspace_transferlearning        - This subproject shows the latent space of chemistry, given that it 
                                                  learns chemical logic, can be used to efficiently transfer learn and solve new problems in chemistry, such as modelling pKa, NMR, and the electron density! We show that pretrained latent spaces of GNN models can be used as inputs to new learning algorithms that are super efficient and require less chemical data. This is essential given that, for some problems, highly reliable chemical data is rare. This also is a test of completeness of GNN models, can they model ALL of chemistry, similar to how the wave function does in quantum chemistry? 


			      

