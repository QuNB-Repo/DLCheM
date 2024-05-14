# DLCheM - What is the chemical logic behing graph neural network (GNN) chemistry


1) gcnn_latentspace_chemistry     - This subproject shows how the latent space of a trained model
                                    is organized on a framework of molecular environments, a.k.a functional groups. 
                                    This organization is so precise, it quantifies the difference between local functional groups with Euclidean distance
                                    This allows the algorithm to finely identify the local structures of chemistry, 
                                    thus basing the decision on this fundamental concept. 
                                    We explore the chemical organization of this space and its ability to quantify molecular similarity


We analyzed two different graph neural networks in chemsitry to arrive at the same result! (hinting that all graph networks of chemistry are doing the same thing!)

    AIMNet_latentspace            - extracts the latent space of a pretrained AIMNet model
    SchNet_latentspace            - extracts the latent space of a pretrained SchNet model
    analysis_tools                - a set of analysis tools to analyze embedddings from either model
                                    dimension reduciton, euclidean distance, classification methods... etc

Please see README within each directory to find out more! 

			      

