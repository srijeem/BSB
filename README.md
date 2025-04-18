# BSB

Pairwise Sequence Alignment using Dynamic Programming (Dynamic_Programming.py):
This Python script implements pairwise sequence alignment using dynamic programming techniques. 
It supports three types of alignments: global, semiglobal, and local. The script can process sequences from a FASTA file and allows users to specify a substitution matrix and gap penalty for the alignment.


Hidden_Markov_Model (Hidden_Markov_Model.py):
The script is a template for implementing Hidden Markov Models (HMMs) in the context of a sequence analysis course. 
It provides basic functionality for training and inference with HMMs, focusing on several core algorithms: 
Viterbi, Forward, Backward, and Baum-Welch. Here's a breakdown of the script's components:

Viterbi Algorithm: This function computes the most probable state sequence for a given observation sequence using dynamic programming. 
It calculates the probability of the sequence and traces back to find the most probable state path.

Forward Algorithm: Similar to the Viterbi algorithm but instead of finding the most probable state sequence, 
it calculates the probability of observing a sequence given the model (forward probability). 
It also uses dynamic programming to build a trellis for intermediate results.

Backward Algorithm: This is the reverse of the Forward algorithm, and is used to calculate the probability of the sequence given the model in a backward manner. 
It computes probabilities recursively starting from the end of the sequence.

Baum-Welch Algorithm: This is an expectation-maximization algorithm used to estimate the parameters of an HMM (transition and emission probabilities) 
based on observed data. 
It iterates between calculating the forward and backward probabilities and updating the parameters to maximize the likelihood of the data.

Main Function: The script's main function loads sequences, transition, and emission matrices, and executes the chosen algorithm
(Viterbi, Forward, Backward, or Baum-Welch) based on user input. 
It handles the input/output of sequences and results and saves relevant output in specified directories.

The script is designed to be run with different commands (viterbi, forward, backward, baumwelch), 
and it can be configured with various verbosity levels to control the amount of output generated during execution. 
Additionally, it provides mechanisms to save the output in specified files, making it easy to analyze the results later.

Score_network_alignments (score_network_alignments.py):

This Python script computes a similarity score for protein network alignments based on Gene Ontology (GO) terms. 
It processes protein-to-Ensembl ID mappings and GO term annotations, calculates the Jaccard index for each protein pair, 
and generates intermediate output files for further analysis. The script is designed to handle protein ID mappings, 
GO term associations, and compute a score quantifying the similarity between two protein networks based on shared GO terms.

Phi-Psi Angles

This Python script processes DSSP (Dihedral Angle and Solvent Accessibility) files to analyze the propensity of different amino acids 
to be buried in the protein core. It computes the fraction of buried surface area (BSA) for each amino acid type and calculates the 
relative propensities for each amino acid to be buried. This information can be useful in the context of protein structure prediction and stability studies.







