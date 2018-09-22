# motifFinder-Bioinformatics

This project involves developing a "motif finding" program and testing the program on a set of synthetic datasets. </br></br> 
The first step is to build a benchmark of the collection of synthetic datasets, and each dataset contains a set of DNA sequences, into which a "motif" has been "planted". </br></br>
The second step is to write a program to read the dataset generated in the first step and find the "motif" planted in the first step. We employed the Gibbs sampling algorithm for implementing the motif finder. </br></br>
Lastly, we evaluated the performance of the motif finder using three metrics: (1) Relative Entropy: difference between the motif our algorithm found and the planted motif, (2) Overlapping Sites: the correct sites our algorithm found, (3) Running time: how long our algorithm takes to find the motif. </br> 
