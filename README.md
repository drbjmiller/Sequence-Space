# Sequence-Space
The Python program perc_sim.py models protein sequence space using one of three possible simulations determined by the value of the variable SIMULATION_TYPE. Each simulation generates a matrix, assigning each cell a random value between 0 and 1. Each sequence corresponds to a cell in the lattice. Cells corresponding to functional sequences hold values lower than the variable Proportion. Each simulation designates a start sequence composed entirely of the first amino acid. The three simulations are the following:
1. Cluster: The "Cluster" simulation recursively follows every path of neighboring functional cells that begins from the star sequence. The simulation records the number of sequences encompassing the entire cluster of paths.
2. Percent: The "Percent" simulation determines if the size of the cluster that includes the start sequence is larger than the variable CLUSTER_MAX. Clusters larger than that value extend throughout the entire sequence space. Nearly all clusters in the trials I ran were smaller than 300 or larger than 500,000. 
3. Attempts: The "Attempts" simulation generates matrices until a continuous path of functional sequences leads from the start sequence to the target defined by a target sequence consisting entirely of the second amino acid and the variable Tol. The simulation records the number of attempts and the path length leading from the start sequence to the target. 

Every simulation accesses the CSV file set by the variable TRIALS_PARAMS_FILE that contains the parameters for simulation trails. Each trial is designed by a row in the file. The parameters include the following: 
* Length: Length of amino acid chain. The length determines the dimensions of a generated matrix that represents sequence space.
* AA Num: Number of amino acids that could reside in each position of the chain. Each dimension has a size of AA Num. Any amino acid can transition to any other amino acid in a single step.
* Proportion: Proportion of sequences that are randomly assigned as functional. 
* Tol: Number of amino acids in a sequence that can differ from a target sequence where the sequence is considered inside the target.
* Steps: The maximum number of amino-acid differences between two sequences that are still considered neighbors.
* Repeat: The number of times a simulation is repeated with the same parameters. 
The simulation launches multiple processes that run in parallel. Each runs all the trials included in TRIALS_PARAMS_FILE. The number of parallel processes is set with the variable PARALLEL_PROC. The number of total trials for each set of parameters is PARALLEL_PROC * Repeat. 

The files trails_params10.csv and trails_params13.csv contained the trails I ran for the matrics corresponding to sequences of length 10 (AA Num = 7) and length 13 (AA Num = 5). The results for those runs were recorded in the file Simulation_Data.xlsx. The simulation adds new results at the beginning of the output files.  

The program only needs to be run with python3. The names of the input and output files are set using the global variables at the beginning of the program. The number of parallel processes that a computer can manage depends on the memory and number of CPU cores. On a Linux server with 512GB RAM and 32 CPU cores, I ran 20 parallel processes with matrices of length 10. Between 5 and 6 parallel processes could be run for matrices of length 13 since they correspond to a much larger sequence space. 
