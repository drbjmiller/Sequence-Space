# Sequence-Space
The Python program perc_sim.py simulates protein sequence space. It accesses the file trials_params.csv which contains the parameters for multiple simulation trails. The parameters include the following: 
* Length: Length of amino acid chain. The length determines the dimensions of a generated matrix that represents sequence space.
* AA Num: Number of amino acids that could reside in each position of the chain. Each dimension has a width of AA Num. Any amino acid can transition to any other amino acid in a single step.
* Proportion: Proportion of sequences that are randomly assigned as functional. 
* Tol: Number of amino acids in a sequence that can differ from a target sequence where the sequence is considered inside the target.
* Steps: The maximum number of amino-acid differences between two sequences that are still considered neighbors.
* Repeat: The number of times a simulation is repeated with the same parameters. 
The simulation launches multiple processes that run in parallel. Each runs all the trials included in trials_params.csv. The number of parallel processes is set with the variable PARALLEL_PROC.

The simulations generate a matrix and assign each cell a random value between 0 and 1. Each sequence corresponds to a cell in the lattice. Cells corresponding to functional sequences hold values lower than the variable Proportion. Each simulation designates a start sequence composed entirely of the first amino acid. The program can run one of three simulations depending on the value of the variable SIMULATION_TYPE. It must be assigned one of three words, which designate the simulation type:
1. Cluster: The "Cluster" simulation assigns a target that determines the target size. The simulation recursively follows every path that begins from the star sequence. It records the number of sequences that encompass the entire cluster of paths.
2. Percent: The "Percent" simulation determines if the size of the cluster that includes the start sequence is larger than the variable CLUSTER_MAX. Clusters larger than that value extend throughout the entire sequence space.
3. Attempts: The "Attempts" simulation generates matrices until a continuous path of functional sequences leads from the start sequence to the target defined by a target sequence consisting entirely of the second amino acid and the variable Tol. The simulation records the number of attempts and the length of the path leading from the start sequence to the target.

The files trails_params10.csv and trails_params13.csv contained the trails I ran for the matrics corresponding to sequences of length 10 (AA Num = 7) and length 13 (AA Num = 5). The results for those runs were recorded in the file Simulation_Data.xlsx. 
