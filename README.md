# Drosophila Gene Clustering Analysis

## Overview

This repository contains the code and documentation for a project focused on clustering patterns of two key genes—cytochrome oxidase subunit I (COI) and alcohol dehydrogenase (Adh)—within the Drosophila genus. The project utilizes unsupervised machine learning techniques, including t-SNE, dendrograms, and silhouette analysis, to explore whether these genes exhibit similar or different clustering behaviors, shedding light on their evolutionary relationships and functional constraints.

## Objectives

- Analyze and compare the clustering patterns of mitochondrial COI and somatic Adh genes within the Drosophila genus.
- Investigate the evolutionary pressures and functional constraints that may influence these clustering patterns.
- Provide insights into the suitability of these genes as molecular markers for classification and evolutionary studies.

## Key Technologies

- **R**: Used for data preprocessing, analysis, and visualization.
- **Bioconductor**: R packages like `msa`, `ape`, and `cluster` for sequence alignment and clustering.
- **t-SNE**: A dimensionality reduction technique used to visualize high-dimensional genetic data.
- **Silhouette Analysis**: A method to measure the quality of clustering.
- **NCBI**: Data source for Drosophila gene sequences.

## Data Description

The sequences for the COI and Adh genes were obtained from the NCBI database, focusing on the Drosophila genus. The study included gene length restrictions and filters to ensure the quality of the sequences, which were then analyzed to identify clustering patterns that could reveal evolutionary and functional relationships within the genus.

## Analysis Workflow

1. **Data Acquisition and Preprocessing**:
   - Retrieved COI and Adh gene sequences from NCBI.
   - Applied gene length restrictions and removed sequences with missing or invalid data.
   - Aligned sequences and removed gaps for cleaner analysis.

2. **Clustering Analysis**:
   - Performed sequence alignment and generated pairwise distance matrices for both genes.
   - Used t-SNE and dendrograms to visualize clustering patterns.
   - Created elbow plots to determine the optimal number of clusters.
   - Conducted silhouette analysis to evaluate clustering strength.

3. **Visualization**:
   - Generated t-SNE plots to compare clustering patterns of COI and Adh genes.
   - Produced silhouette plots to assess the quality of clustering for each gene.

## Results

- **Clustering Patterns**: The t-SNE analysis revealed that COI and Adh genes cluster differently within the Drosophila genus, indicating diverse evolutionary paths.
- **Elbow and Silhouette Analysis**: The optimal number of clusters for COI was determined to be 5, while for Adh it was 3. The silhouette analysis showed stronger clustering for Adh compared to COI, with some misclassified clusters in the COI gene.
- **Evolutionary Insights**: The differences in clustering patterns suggest that COI and Adh are subject to distinct evolutionary pressures and functional roles, with implications for their use as genetic markers in evolutionary studies.

## Future Directions

- Expand the analysis to include additional genes from the Drosophila genus to validate findings.
- Investigate the impact of using mitochondrial versus somatic genes for clustering analysis.
- Explore other machine learning techniques to improve clustering accuracy and interpretability.

## References

- Peng, H., Long, F., Eisen, M. B., & Myers, E. W. (2006). Clustering gene expression patterns of fly embryos. *IEEE International Symposium on Biomedical Imaging*.
- Leger, J., Daudin, J., & Vacher, C. (2015). Clustering methods differ in their ability to detect patterns in ecological networks. *Methods in Ecology and Evolution*.
- Liu, C., Wright, B., Allen-Vercoe, E., Gu, H., & Beiko, R. (2018). Phylogenetic clustering of genes reveals shared evolutionary trajectories and putative gene functions. *Genome Biology and Evolution*.

