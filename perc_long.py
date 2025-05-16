#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 15:52:21 2022
Updated on Monday Nov 27 2023

@author: Brian Miller

This program simulates protein sequence space by exploring an array whose dimension equals the length
of a protein sequence. Each dimension has a length equal to the number of amino acids: 20. Amino acids 
can transition to other amino acids based on the genetic code. Unlike the original version of the code,
the entire sequence space is not created in advance, but sampled sequences are stored in a list. For 
each percentage of functional sequences, the simulation is run 100 times. The simulation records the 
size of the cluster containing the initial sequence, and it records the percentage of clusters larger
than a threshold for each set of 100 trials. 

"""
import pandas as pd
import random
import sys

COLUMNS_PERCENT = ['Length', 'Pfs', 'Offset', 'Percent'] 
COLUMNS_CLUSTERS = ['Length', 'Pfs', 'Offset', 'Size'] 

CLUSTER_THRESHOLD = 2000
AA_TRANS = [[8, 11, 15, 16, 17, 18],
            [4, 13, 17, 18, 19],
            [3, 7, 8, 10, 11, 16, 19],
            [2, 6, 7, 10, 12, 13, 16],
            [1, 9, 10, 13, 17, 19],
            [6, 7, 8, 14, 17, 18],
            [3, 5, 8, 12, 13, 16],
            [2, 3, 5, 10, 14, 17, 18],
            [0, 2, 5, 6, 11, 15, 17],
            [4, 10, 11, 16, 18, 19],
            [2, 3, 4, 7, 9, 19],
            [0, 2, 8, 9, 15, 16, 17, 18, 19],
            [3, 6, 13, 14, 15, 16, 19],
            [1, 3, 4, 6, 12, 16, 17, 19],
            [5, 7, 12, 15, 17, 18, 19],
            [0, 2, 8, 11, 12, 14, 17, 19],
            [1, 3, 6, 9, 11, 12, 13, 18],
            [0, 1, 4, 5, 7, 8, 11, 13, 14, 15, 18, 19],
            [0, 1, 5, 7, 9, 11, 14, 16, 17],
            [1, 2, 4, 9, 19, 11, 12, 13, 14, 15, 17, 18]]

LENGTHS = [500]

# Iterative version of the cluster size determination
def cluster_size(length, start_seq, Pfs):
    traversed_seqs = set()
    queue = [start_seq]
    cluster_size = 0

    while queue:
        if cluster_size > CLUSTER_THRESHOLD:
            break
        current_seq = queue.pop()
        traversed_seqs.add(tuple(current_seq))
        cluster_size += 1
        if cluster_size % 500 == 0:
            print("Cluster Size: {}".format(cluster_size))

        for pos_change in range(length):
            possible_aa = AA_TRANS[current_seq[pos_change]]
            for new_aa in possible_aa:
                new_seq = current_seq.copy()
                new_seq[pos_change] = new_aa
                func_check = random.random()

                if func_check > Pfs:
                    traversed_seqs.add(tuple(new_seq))
                    continue

                if tuple(new_seq) not in traversed_seqs:
                    queue.append(new_seq)

    return cluster_size

# Main Program
sys.setrecursionlimit(15000)
data_percent = pd.DataFrame(columns=COLUMNS_PERCENT)
data_clusters = pd.DataFrame(columns=COLUMNS_CLUSTERS)

START_OFFSET = 0
END_OFFSET = .2

for length in LENGTHS:
    Pth = 1 / (7.5 * length)
    offset = START_OFFSET
    large_percent = 0

    while offset < END_OFFSET - .01:
        Pfs = Pth * (1 + offset)
        large_clusters = 0

        for trial in range(100):
            if trial % 5 == 0 and trial != 0:
                print("Trial: {}".format(trial))
            initial_seq = [random.randint(0, 19) for _ in range(length)]
            size = cluster_size(length, initial_seq, Pfs)
            data_clusters.loc[len(data_clusters.index)] = [length, Pfs, offset, size]

            if size >= CLUSTER_THRESHOLD:
                large_clusters += 1

        large_percent = large_clusters / 100
        data_percent.loc[len(data_percent.index)] = [length, Pfs, offset, large_percent]
        print("Length: {}  Offset: {:.2f}  Pfs: {:.6f}  Large Percent: {:.2f}".format(length, offset, Pfs, large_percent))
        offset += 0.02
        
# Outputs data to files
file_percent = "data_percent " + str(LENGTHS[0]) + " " + str(END_OFFSET) + ".csv"
file_size = "data_size " + str(LENGTHS[0]) + " " + str(END_OFFSET) + ".csv"
data_percent.to_csv(file_percent, encoding='utf-8', index=False) 
data_clusters.to_csv(file_size, encoding='utf-8', index=False) 

