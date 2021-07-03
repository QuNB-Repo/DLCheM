# DLCheM

Before starting, please review the [contribution guidelines](CONTRIBUTING.md). The (Deep) Machine Learning Chemistry repo aims to answer the following question: [Can we learn chemistry by interpreting the Schnet deep neural network](DOCS.md)? 

## Directories

schnet_project --> directory for work in investigating the interpretability of schnet algorithm using statistical techniques such as principle component analysis (PCA). Within this directory you will find: 

	notebooks - folder of jupyter files for running schnet training script, load trained model, and sklearn PCA   
	            script.
	pca -       folder containing all PCA results for the many trained models 
	schnet-package-install 
	trained_models - folder which saves all the trained_models after training is done

bpy-battery_project --> work in searching for optimal bpy battery properties using machine learning algorithms. 			within this folder you will find: 

	autocombinatorial - folder containing the script for autocombinatorial generation of new bpy molecules 
	SMILES-RNN -        working on understanding SMILES representations in ML algorithms

