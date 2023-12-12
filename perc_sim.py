#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 15:52:21 2022
Updated on Monday Nov 27 2023

@author: Brian Miller

This program simulates protein sequence space by creating matrices whose dimension equals the length
of a peptide sequence. Each dimension has a length equal to the number of amino acids. The value 
of each site is set to a random number between 0 and 1. If the value is smaller than a set proportion of 
functional sequences, the sequence corresponding to that site is considered functional. 

The simulation determines the following:
    * If SIMULATION_TYPE = "Cluster", the simulation reports the size of the cluster containing the start sequence.
    * If SIMULATION_TYPE = "Percent", the simulation reports the percentage of start sequences in large clusters
      extending throughout sequence space.
    * If SIMULATION_TYPE = "Attempts", the simulation reports the number of attempts 
      required for a path to connect the start sequence to the target.
The variable "steps" represents the maximum number of amino acid differences between two sequences that are still neighbors. 
It can equal 1, 2, or 3. 
"""
import multiprocessing
import numpy as np
import pandas as pd
import random
import sys
import os.path
from os import getpid
import math 

COLUMNS_CLUSTER = ['Length', 'AA Num', 'Process', 'PID', 'Proportion', 'Cluster'] 
COLUMNS_CLAV = ['Length', 'AA Num', 'Proportion', 'Cluster Ave', 'STD', 'Trials'] 
COLUMNS_PERC_LARGE = ['Length', 'AA Num', 'Process', 'PID', 'Proportion', 'Size'] 
COLUMNS_PLGAV = ['Length', 'AA Num', 'Proportion', 'Perc Large', 'STD', 'Trials'] 
COLUMNS_ATTEMPTS = ['Length', 'AA Num', 'Process', 'PID', 'Proportion', 'Attempts', 'Path Len'] 
COLUMNS_ATTAV = ['Length', 'AA Num', 'Proportion', 'Attempts Ave', 'Att STD', 'Path Len Ave', 'Len STD', 'Trials'] 

PARALLEL_PROC = 2                               # Number of parallel processes
PATH_MAX = 15000                                # Maximum path length for average attempts and average cluster size
PATH_MAX_PLG = 3000                             # Maximum attempts for % large clusters
CLUSTER_MAX = 20000                             # Maximum cluster size for % large clusters

SIMULATION_TYPE = "Cluster"                     # Type of simulation: "Cluster", "Percent", or "Attempts"
TRIALS_PARAMS_FILE = "trials_params.csv"        # Name of file with parameter values for trials
PROCESS_FILE_BASE = "process"                   # File name base for output of processes
CL_OUTPUT_FILE = "output_cl.csv"                # Name of output file for individual cluster output
CLAV_OUTPUT_FILE = "output_clav.csv"            # Name of output file for average cluster output
PLG_OUTPUT_FILE = "output_plg.csv"              # Name of output file for individual % large cluster output
PLGAV_OUTPUT_FILE = "output_plgav.csv"          # Name of output file for average % large cluster output
ATT_OUTPUT_FILE = "output_att.csv"              # Name of output file for individual number of attempts
ATTAV_OUTPUT_FILE = "output_attav.csv"          # Name of output file for average number of attempts


#############################################
# Create new data file if needed
def test_datafile(file_name, column_names):
    df = pd.DataFrame(columns = column_names)
    # Check if file exists. If not, create it.
    if os.path.exists(file_name) == False:
        df.to_csv(file_name, encoding='utf-8', index=False)
        return
    
    # check if file column names are correct. If not, create new file with correct names.
    data = pd.read_csv(file_name)
    if data.columns.values.tolist() != column_names:
        df.to_csv(file_name, encoding='utf-8', index=False)
        
    return


#############################################
# Create grid of dimensions length and size aa with random values
def initialize(length, aa):
    if length == 7:
        mat = np.random.rand(aa,aa,aa,aa,aa,aa,aa)
        mat[1][1][1][1][1][1][1] = 0
    elif length == 8:
        mat = np.random.rand(aa,aa,aa,aa,aa,aa,aa,aa)
        mat[1][1][1][1][1][1][1][1] = 0
    elif length == 9:
        mat = np.random.rand(aa,aa,aa,aa,aa,aa,aa,aa,aa)
        mat[1][1][1][1][1][1][1][1][1] = 0
    elif length == 10:
        mat = np.random.rand(aa,aa,aa,aa,aa,aa,aa,aa,aa,aa)
        mat[1][1][1][1][1][1][1][1][1][1] = 0
    elif length == 11:
        mat = np.random.rand(aa,aa,aa,aa,aa,aa,aa,aa,aa,aa,aa)
        mat[1][1][1][1][1][1][1][1][1][1][1] = 0
    elif length == 12:
        mat = np.random.rand(aa,aa,aa,aa,aa,aa,aa,aa,aa,aa,aa,aa)
        mat[1][1][1][1][1][1][1][1][1][1][1][1] = 0
    elif length == 13:
        mat = np.random.rand(aa,aa,aa,aa,aa,aa,aa,aa,aa,aa,aa,aa,aa)
        mat[1][1][1][1][1][1][1][1][1][1][1][1][1] = 0
    else:
        print("Dimension %d not included" % length)
        sys.exit("\n")
    return(mat)


#############################################
# Check if sequences differ by no more than tol
def check_same(list1, list2, tol):
    same = True
    dif = 0
    for i in range(len(list1)):
        if list1[i] != list2[i]:
            dif += 1
    if dif > tol:
        same = False
    return(same)
        

#############################################
# Retrieve value in matrix at position corresponding to seq
def retrieve(mat, seq, length):
    if length == 1:
        value = mat[seq[0]] 
        return value
    else:
        return(retrieve(mat[seq[0]], seq[1:], length-1))


#############################################
# Set value of cell designated by seq to 1 to make cell impassable
def block_path(mat, seq, length):
    if length == 1:
        mat[seq[0]] = 1
    else:
        block_path(mat[seq[0]], seq[1:], length-1)
    
    
#############################################
# Recursively traverse paths through grid
def search_path(mat, length, aanum, proportion, steps, tol, seq, target_seq, pathlen, pathinfo):
    # Block sequence from becoming part of other paths. 
    block_path(mat, seq, length)
    pathlen += 1
    
    # Check if path length > PATH_MAX. Simulation crashes if path becomes too large. 
    if pathlen > PATH_MAX:
        return
    
    # Check type of simulation to determine next actions
    if SIMULATION_TYPE == "Cluster":
        pathinfo["Cluster"] += 1        # Increment cluster size by 1
    elif SIMULATION_TYPE == "Percent":
        # Check if cluster size exceeds CLUSTER_MAX indicating large cluster.
        if pathinfo["Size"] == "Large":
            return
        pathinfo["Cluster"] += 1        # Increment cluster size by 1
        # Check if cluster size exceeds CLUSTER_MAX. If so, set pathinfo["Cluster Max"] to True
        if pathinfo["Cluster"] >= CLUSTER_MAX:
            pathinfo["Size"] = "Large"
            return
    elif SIMULATION_TYPE == "Attempts":
        if pathinfo["Path Found"] == True:
            return
        # Check it target found. If found, record pathlen and end all searches
        if check_same(seq, target_seq, tol):
            pathinfo["Path Found"] = True
            pathinfo["Path Len Target"] = pathlen
            return
    
    
    # Explore different paths to see if a functional sequence resides within "steps" mutations
    aas = [item for item in range(aanum)]
    if steps == 1:
        for pos_change in range(length):
            for newaa in random.sample(aas, aanum):
                newseq = seq.copy()
                newseq[pos_change] = newaa
                if retrieve(mat, newseq, len(newseq)) < proportion:
                    search_path(mat, length, aanum, proportion, steps, tol, newseq, target_seq, pathlen, pathinfo)
    elif steps == 2:
        for pos_change1 in range(length):
            for newaa1 in random.sample(aas, aanum):
                for pos_change2 in range(pos_change1, length):
                    for newaa2 in random.sample(aas, aanum):
                        newseq = seq.copy()
                        newseq[pos_change1] = newaa1
                        newseq[pos_change2] = newaa2                   
                        if retrieve(mat, newseq, len(newseq)) < proportion:
                            if pos_change1 != pos_change2:
                                search_path(mat, length, aanum, proportion, steps, tol, newseq, target_seq, pathlen+2, pathinfo)
                            else:
                                search_path(mat, length, aanum, proportion, steps, tol, newseq, target_seq, pathlen+1, pathinfo)
    elif steps == 3:
        for pos_change1 in range(length):
            for newaa1 in random.sample(aas, aanum):
                for pos_change2 in range(pos_change1, length):
                    for newaa2 in random.sample(aas, aanum):
                        for pos_change3 in range(pos_change2, length):
                            for newaa3 in random.sample(aas, aanum):
                                newseq = seq.copy()
                                newseq[pos_change1] = newaa1
                                newseq[pos_change2] = newaa2                   
                                newseq[pos_change3] = newaa3                   
                                if retrieve(mat, newseq, len(newseq)) < proportion:
                                    if pos_change1 != pos_change2 and pos_change2 != pos_change3:
                                        search_path(mat, length, aanum, proportion, steps, tol, newseq, target_seq, pathlen+3, pathinfo, True)
                                    elif pos_change1 != pos_change2 or pos_change2 != pos_change3:
                                        search_path(mat, length, aanum, proportion, steps, tol, newseq, target_seq, pathlen+2, pathinfo, True)
                                    else:
                                        search_path(mat, length, aanum, proportion, steps, tol, newseq, target_seq, pathlen+1, pathinfo, True)
    else:
        print("Steps: %d not included" % steps)
        sys.exit("\n")

    return

#############################################
# Run separate processes to search all p values PARALLEL_PROC times in parallel
def perc_process(process_num):
    if SIMULATION_TYPE == "Cluster":
        process_file = PROCESS_FILE_BASE + "_cl" + str(process_num) + ".csv"
        test_datafile(process_file, COLUMNS_CLUSTER)
        pathinfo = {"Cluster": 0}
    elif SIMULATION_TYPE == "Percent":
        process_file = PROCESS_FILE_BASE + "_plg" + str(process_num) + ".csv"
        test_datafile(process_file, COLUMNS_PERC_LARGE)
        pathinfo = {"Cluster": 0, "Size": "Small"}
    elif SIMULATION_TYPE == "Attempts":
        process_file = PROCESS_FILE_BASE + "_att" + str(process_num) + ".csv"
        test_datafile(process_file, COLUMNS_ATTEMPTS)
        pathinfo = {"Path Found": False, "Path Len Target": 0}        
    
    pid = getpid()
    
    # Iterate through trials in trials_params.csv
    trials_params = pd.read_csv(TRIALS_PARAMS_FILE)
    for index, row in trials_params.iterrows():
        # Read run parameters from file
        length = int(row['Length'])
        aanum = int(row['AA Num'])
        proportion = float(row['Proportion'])
        tol = int(row['Tol'])
        steps = int(row['Steps'])
        repeat = int(row['Repeat'])

        # Create initial and target sequences
        initial_seq = [0]*length
        target_seq = [1]*length
        
        # Run traials for same paramaeters repeat times.
        for run in range(repeat):
            if SIMULATION_TYPE == "Cluster":
                # Initialize variables and matrix
                pathlen = -1
                pathinfo["Cluster"] = 0
                mat = initialize(length, aanum)
                
                # Search through matricies for all paths
                search_path(mat, length, aanum, proportion, steps, tol, initial_seq, target_seq, pathlen, pathinfo)
        
                # Record results to file
                trial_data = [length, aanum, process_num, pid, proportion, pathinfo["Cluster"]]
                process_data = pd.read_csv(process_file)
                process_data.loc[len(process_data.index)] = trial_data
                process_data[['Length', 'AA Num', 'Process', 'PID', 'Cluster']] = process_data[['Length', 'AA Num', 'Process', 'PID', 'Cluster']].astype('int64')
                process_data.to_csv(process_file, encoding='utf-8', index=False)
            elif SIMULATION_TYPE == "Percent":
                # Initialize variables and matrix
                pathlen = -1
                pathinfo["Cluster"] = 0
                pathinfo["Size"] = "Small"                
                mat = initialize(length, aanum)
                
                # Search through matricies for all paths
                search_path(mat, length, aanum, proportion, steps, tol, initial_seq, target_seq, pathlen, pathinfo)
    
                # Record results to file
                process_data = [length, aanum, process_num, pid, proportion, pathinfo["Size"]]
                perc_data = pd.read_csv(process_file)
                perc_data.loc[len(perc_data.index)] = process_data
                perc_data[['Length', 'AA Num', 'Process', 'PID']] = perc_data[['Length', 'AA Num', 'Process', 'PID']].astype('int64')
                perc_data.to_csv(process_file, encoding='utf-8', index=False)
                
            elif SIMULATION_TYPE == "Attempts":
                # Initialize variables
                pathinfo["Path Found"] = False
                pathinfo["Path Len Target"] = 0
                attempts = 0

                # Search through matricies until target found 
                while pathinfo["Path Found"] == False:
                    pathlen = -1                      # initialize path length
                    mat = initialize(length, aanum)   # Initialize matrix
                    search_path(mat, length, aanum, proportion, steps, tol, initial_seq, target_seq, pathlen, pathinfo)
                    attempts += 1
        
                # Record results to file
                output_data = [length, aanum, process_num, pid, proportion, attempts, pathinfo["Path Len Target"]]
                perc_data = pd.read_csv(process_file)
                perc_data.loc[len(perc_data.index)] = output_data
                perc_data[['Length', 'AA Num', 'Process', 'PID', 'Attempts', 'Path Len']] = perc_data[['Length', 'AA Num', 'Process', 'PID', 'Attempts', 'Path Len']].astype('int64')
                perc_data.to_csv(process_file, encoding='utf-8', index=False)                
    

################### Main Program ##########################

sys.setrecursionlimit(PATH_MAX + 1000)
simulation_types = ['Cluster', 'Percent', 'Attempts']
if SIMULATION_TYPE not in simulation_types:
    sys.exit("You did not choose a valid simulation type. \nSIMULATION_TYPE must equal \"Cluster\", \"Percent\", or \"Attempts\"")

if __name__ == '__main__':
    processes = []
    # Create PARALLEL_PROC processes for parallel computing
    for process_num in range(PARALLEL_PROC):
        process = multiprocessing.Process(target=perc_process, args=[process_num+1])
        processes.append(process)
        process.start()
        
    # Combine processes files into single output file. 
    process_num = 0
    frames = []
    for process in processes:
        process_num += 1
        process.join()
        if SIMULATION_TYPE == "Cluster":
            proc_data_file = PROCESS_FILE_BASE + "_cl" + str(process_num) + ".csv"
        elif SIMULATION_TYPE == "Percent":
            proc_data_file = PROCESS_FILE_BASE + "_plg" + str(process_num) + ".csv"
        elif SIMULATION_TYPE == "Attempts":
            proc_data_file = PROCESS_FILE_BASE + "_att" + str(process_num) + ".csv"
        perc_data = pd.read_csv(proc_data_file)
        frames.append(perc_data)
        os.remove(proc_data_file)

    # Save data to the output file and open/create a file to write the average cluster size, percentage of large clusters, or average attempts
    if SIMULATION_TYPE == "Cluster":
        test_datafile(CL_OUTPUT_FILE, COLUMNS_CLUSTER)
        output_old_data = pd.read_csv(CL_OUTPUT_FILE)
        frames.append(output_old_data)
        output_data = pd.concat(frames)
        output_data.to_csv(CL_OUTPUT_FILE, encoding='utf-8', index=False)
    elif SIMULATION_TYPE == "Percent":
        test_datafile(PLG_OUTPUT_FILE, COLUMNS_PERC_LARGE)
        output_old_data = pd.read_csv(PLG_OUTPUT_FILE)
        frames.append(output_old_data)
        output_data = pd.concat(frames)
        output_data.to_csv(PLG_OUTPUT_FILE, encoding='utf-8', index=False)
    elif SIMULATION_TYPE == "Attempts":
        test_datafile(ATT_OUTPUT_FILE, COLUMNS_ATTEMPTS)
        output_old_data = pd.read_csv(ATT_OUTPUT_FILE)
        frames.append(output_old_data)
        output_data = pd.concat(frames)
        output_data.to_csv(ATT_OUTPUT_FILE, encoding='utf-8', index=False) 

    # Calculate average cluster size, percent large clusters, or average attempts and same to output file. 
    if SIMULATION_TYPE == "Cluster":
        output_ave = pd.DataFrame(columns = COLUMNS_CLAV)     
    elif SIMULATION_TYPE == "Percent":
        output_ave = pd.DataFrame(columns = COLUMNS_PLGAV) 
    elif SIMULATION_TYPE == "Attempts":        
        output_ave = pd.DataFrame(columns = COLUMNS_ATTAV) 

    length_set = list(set(output_data.Length))                      # List of sequence lengths in output data
    for length in length_set:
        output_length = output_data[output_data.Length==length]
        proportion_set = list(set(output_length.Proportion))        # List of proportions for each length
        for proportion in proportion_set:
            output_len_prop = output_length[output_length.Proportion==proportion]
            sample = len(output_len_prop)
            aa_num = output_len_prop.iloc[0]['AA Num']
            # Average relevant variable(s) for simulation type
            if SIMULATION_TYPE == "Cluster":
                cluster_av = round(output_len_prop["Cluster"].mean(),1)      # Calculate average cluster size
                cluster_std = round(output_len_prop["Cluster"].std()/math.sqrt(sample),1)     # Calculate average cluster size
                output_avdata = [length, aa_num, proportion, cluster_av, cluster_std, sample]
                output_file = CLAV_OUTPUT_FILE
            elif SIMULATION_TYPE == "Percent":
                # Calculate percentage of large clusters
                percent_lg = round(len(output_len_prop[output_len_prop.Size == "Large"]) / sample, 2)
                perclg_std = round(math.sqrt(percent_lg * (1 - percent_lg) / sample), 3)
                output_avdata = [length, aa_num, proportion, percent_lg, perclg_std, sample]
                output_file = PLGAV_OUTPUT_FILE
            elif SIMULATION_TYPE == "Attempts":
                attempts_ave = round(output_len_prop["Attempts"].mean(), 1)
                attempts_std = round(output_len_prop["Attempts"].std()/math.sqrt(sample), 1)
                path_len_ave = round(output_len_prop["Path Len"].mean(), 1) 
                path_len_std = round(output_len_prop["Path Len"].std()/math.sqrt(sample) , 1)
                output_avdata = [length, aa_num, proportion, attempts_ave, attempts_std, path_len_ave, path_len_std, sample]
                output_file = ATTAV_OUTPUT_FILE
            output_ave.loc[len(output_ave.index)] = output_avdata   
            output_ave[['Length', 'AA Num', 'Trials']] = output_ave[['Length', 'AA Num', 'Trials']].astype('int64')      
    output_ave.to_csv(output_file, encoding='utf-8', index=False)     

    print("All percolation processes have finished")

