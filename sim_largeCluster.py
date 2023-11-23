#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 15:52:21 2022

@author: bmiller

Simulation determines if start sequence is part of cluster that extends throughout sequence space.
ns can equal 1, 2, or 3
"""
import multiprocessing
import numpy as np
import pandas as pd
import random
import sys
import os.path
from os import getpid
import csv
from statistics import mean 

COLUMN_NAMES_DATA = ['Length', 'AA Num', 'Process', 'PID', 'Proportion', 'Transition'] 

PARALLEL_PROC = 2
PATH_MAX = 3000
CLUSTER_MAX = 20000
OUTPUT_FILE_BASE = "process"
OUT_FILE_COMBINED = "out_lgCluster.csv"
TRIALS_PARAMS_FILE = "trials_params.csv"


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
def exclude(mat, seq, length):
    if length == 1:
        mat[seq[0]] = 1
    else:
        exclude(mat[seq[0]], seq[1:], length-1)
    
    
#############################################
# Recursively traverse paths through grid
def search_path(mat, length, aanum, proportion, steps, tol, seq, pathlen, pathinfo):
    exclude(mat, seq, length)
    pathlen += 1
    pathinfo["Cluster"] += 1
    
    # Check if path length > PATH_MAX or cluster size exceeded CLUSTER_MAX. If either is true, end path.
    if pathlen >= PATH_MAX or pathinfo["Cluster Max"] == True:
        return

    # Check if cluster size exceeds CLUSTER_MAX. If so, set pathinfo["Cluster Max"] to True
    if pathinfo["Cluster"] >= CLUSTER_MAX:
        pathinfo["Cluster Max"] = True
        return

    # Explore different paths to see if a functional sequence resides within "steps" mutations
    aas = [item for item in range(aanum)]
    if steps == 1:
        for pos_change in range(length):
            for newaa in random.sample(aas, aanum):
                newseq = seq.copy()
                newseq[pos_change] = newaa
                if retrieve(mat, newseq, len(newseq)) < proportion:
                    search_path(mat, length, aanum, proportion, steps, tol, newseq, pathlen, pathinfo)
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
                                search_path(mat, length, aanum, proportion, steps, tol, newseq, pathlen+2, pathinfo)
                            else:
                                search_path(mat, length, aanum, proportion, steps, tol, newseq, pathlen+1, pathinfo)
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
                                        search_path(mat, length, aanum, proportion, steps, tol, newseq, pathlen+3, pathinfo)
                                    elif pos_change1 != pos_change2 or pos_change2 != pos_change3:
                                        search_path(mat, length, aanum, proportion, steps, tol, newseq, pathlen+2, pathinfo)
                                    else:
                                        search_path(mat, length, aanum, proportion, steps, tol, newseq, pathlen+1, pathinfo)
    else:
        print("Steps: %d not included" % steps)
        sys.exit("\n")

    return

#############################################
# Run separate processes to search all p values multiple times in parallel
def perc_process(process_num):
    process_file = OUTPUT_FILE_BASE + str(process_num) + ".txt"
    test_datafile(process_file, COLUMN_NAMES_DATA)
    pathinfo = {"Cluster":0, "Cluster Max": False}
    np.random.seed(process_num)
    pid = getpid()
    
    # Interate through trials in perc_trials.csv
    trials_params = pd.read_csv(TRIALS_PARAMS_FILE)
    for index, row in trials_params.iterrows():
        # Read run parameters from file
        length = int(row['Length'])
        aanum = int(row['AA Num'])
        proportion = float(row['Proportion'])
        tol = int(row['Tol'])
        steps = int(row['Steps'])
        runs = int(row['Runs'])

        # Initialize variables
        initial_seq = [0]*length
        
        for run in range(runs):
            
            # Search through matricies for all paths
            pathlen = -1
            pathinfo["Cluster"] = 0
            pathinfo["Cluster Max"] = False
            mat = initialize(length, aanum)
            search_path(mat, length, aanum, proportion, steps, tol, initial_seq, pathlen, pathinfo)

            # Record results to file
            process_data = [length, aanum, process_num, pid, proportion, pathinfo["Cluster Max"]]
            perc_data = pd.read_csv(process_file)
            perc_data.loc[len(perc_data.index)] = process_data
            perc_data[['Length', 'AA Num', 'Process', 'PID']] = perc_data[['Length', 'AA Num', 'Process', 'PID']].astype('int64')
            perc_data.to_csv(process_file, encoding='utf-8', index=False)
    

################### Main Program ##########################

sys.setrecursionlimit(PATH_MAX + 1000)

if __name__ == '__main__':
    processes = []
    for process_num in range(PARALLEL_PROC):
        process = multiprocessing.Process(target=perc_process, args=[process_num+1])
        processes.append(process)
        process.start()
        
    process_num = 0
    frames = []
    for process in processes:
        process_num += 1
        process.join()
        perc_data_file = OUTPUT_FILE_BASE + str(process_num) + ".txt"
        perc_data = pd.read_csv(perc_data_file)
        frames.append(perc_data)
      
    output_data = pd.concat(frames)
    output_data.to_csv(OUT_FILE_COMBINED, encoding='utf-8', index=False)
    print("All percolation processes have finished")

