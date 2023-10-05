# DLCheM - Can we learn chemistry from interpreting SchNet Neural Networks? 

What sorts of chemistry is hidden inside statistical representation of trained neural networks? 

If we train a SchNet graph convolutional algorithm, why does it make the decisions that it does? 

To answer this question, we extract the nodes in the graph that SchNet builds for each atom (each node corresponding to an atom). These nodes contain the information
that SchNet learned in what we call "embedding vectors." 1) We found these vectors to contain a ton of fundamental chemistry concepts such as a notion of functional group chemistry, and 
a notion of molecular similarity, and much more. 2) We also found that this representation is ideal to transfer learn and predict properties of totally new chemical properties (such as pKa, NMR, 
solubility) using little data and very rudimentary statistical methods.

Inside this directory there are two main projects:

	1) extractembeddings - this code can extract embedding vectors, the representation of atom-in-molecules built by SchNet's graph convolutional neural network.  
                               it also extracts other relevant information with the embedding such as: atom functional group label and colors, interaction layer #, molecule index, element.
			       read the README file in extractembeddings to see how the output file is built.

 
 	2) transferlearning - this code contains simple statistical methods (linear regression, feedforward neural nets, and atomwise size-extensive neural nets) that can transfer learn, 
                             from embedding representation of atoms to atomistic or molecular properties.
			      

