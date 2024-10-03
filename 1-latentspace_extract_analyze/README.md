# Latent Space Extraction of Various Deep Learning Models

This project focuses on extracting the latent space of trained models 
and analyzing its organization based on molecular environments, 
specifically functional groups. By reducing the dimensionality of 
the latent space, we can observe how it captures the structural 
similarities between functional groups, quantified by Euclidean 
distance. This enables the algorithm to identify and classify 
local molecular structures in chemistry, making decisions based 
on these fundamental concepts.

---

### 2) **Subprojects and Tools**

- `AIMNet_latentspace`: Extracts the latent space of a pretrained 
  AIMNet model.
- `SchNet_latentspace`: Extracts the latent space of a pretrained 
  SchNet model.
- `analysis_tools`: A set of analysis tools to evaluate embeddings 
  from either model. These tools include dimension reduction 
  techniques, Euclidean distance calculations, classification 
  methods, and more.

---

## Latent Space Analysis

By extracting and analyzing the latent space of both AIMNet and 
SchNet, this project demonstrates that the organization of 
molecular environments is a fundamental characteristic across 
different graph neural network models in chemistry. The precision 
of these models in identifying functional groups and molecular 
similarities through Euclidean distance highlights the importance 
of latent space analysis in understanding chemical structure.

---

## Tools and Methods

This project utilizes various methods for analyzing latent space, 
including:
- **Dimension reduction**: Techniques like t-SNE and PCA are used 
  to visualize the high-dimensional latent space.
- **Euclidean distance**: A key metric to quantify the similarity 
  between local functional groups.
- **Classification**: Tools to classify functional groups based on 
  their organization in the latent space.

---

For more detailed explanations, refer to the README files within 
each subproject directory.
