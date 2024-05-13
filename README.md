# DLCheM - What is the chemical logic behing graph neural network (GNN) chemistry

How do we analyze "black box" models, specifically grah neural networks (GNN) models on chemistry? Does the information they store after training hold any meaningful trends? Is it possible to break down the decision making of such models, and if so, is it using proper chemical logic as part of its decision-making framework? These are the questions that must be tackled if AI is to revolutionalize the chemical sciences with its speed, efficiency, and trustworthy accurate predictions.

This project analyzes the latent space of graph convolutional neural networks, in particular SchNet and AIMNet, 
to reveal their connection to chemical logic and chemical syntax, and how this logic learned by the algorithm can be used to transfer solve new problems in chemistry very efficiently! 

The project is divided in three subprojects: 

    1) gcnn_latentspace_chemistry               - This subproject shows how the latent space of a trained GNN model
                                                  on chemistry is a space that is organized based on a framework of molecular environments, a.k.a functional groups. It is so precise, it quantifies the difference between local functional groups. This allows the algorithm to identify the local structures of chemistry, thus basing the ultimate decision on this fundamental concept. We explore the chemical organization of this space and its ability to quantify molecular similarity between chemical shapes
    2) gcnn_latentspace_reactions               - This subproject shows how the latent space of a trained model obeys 
                                                  chemical reaction syntax. Similar to how natural language algorithms can solve word analogies, in a famous example "King" - "Man" + "Woman" is nearest in latent space to vector for "Queen". Here we run reactions in a database to find out the same can be said about chemistry GNN models, "amides" - "amine" + "alcohol" is nearest in the latent space to "carboxylic acid" vector. This reveals that the latent space is organized based on reaction syntax, we can drive all reactants of a certain kind to their products with a single, constanct vector, for all the reactants. This is true for any reaction chosen. The algorithm in this project is able to to run any reactions in a database and analyze the embeddings thereof to reveal these chemical reaction analogies in the latent space.
    3) gcnn_latentspace_transferlearning        - This subproject shows the latent space of chemistry, given that it 
                                                  learns chemical logic, can be used to efficiently transfer learn and solve new problems in chemistry, such as modelling pKa, NMR, and the electron density! We show that pretrained latent spaces of GNN models can be used as inputs to new learning algorithms that are super efficient and require less chemical data. This is essential given that, for some problems, highly reliable chemical data is rare. This also is a test of completeness of GNN models, can they model ALL of chemistry, similar to how the wave function does in quantum chemistry? 
    4) flow_batteries_bpy_genetic algorithm     - an application of genetic algorithms on bipyridine derivatives to 
                                                  predict a bipyridine-like substitute that has better flow battery properties as a method of inverse molecular design.

    See README files within each subproject to learn more! 


			      

