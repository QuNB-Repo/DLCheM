# DLCheM - Can we learn chemistry from interpreting SchNet Neural Networks? 

 Inside this directory there are two main projects:

	1) extractembeddings - this code can extract embedding vectors, the representation of atom-in-molecules built by SchNet's graph convolutional neural network.  
                               it also extracts other relevant information with the embedding such as: atom functional group label and colors, interaction layer #, molecule index, element.
			       read the README file in extractembeddings to see how the output file is built.

 
 	2) transferlearning - this code contains simple statistical methods (linear regression, feedforward neural nets, and atomwise size-extensive neural nets) that can transfer learn, 
                             from embedding representation of atoms to atomistic or molecular properties. Read the internal README file for this project for more information
			      

