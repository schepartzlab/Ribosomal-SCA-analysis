# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 18:20:27 2018

@author: Allison Walker
This script calculates the coupling between all pairs of sectors and groups highly coupled sectors from an SCA output file
Cutoff used to group a pair of sectors is .256, which was the mean of the coupling between all pairs of sectors + 2 standard deviations

Prints a list of the grouped sectors
Writes 3 files:
    between_sector_coupling_avg.txt a matrix containing the average coupling between sectors
    between_sector_coupling_max.txt a matrix containing the maximum coupling between sectors
    between_sector_coupling_median.txt a matrix containing the median coupling between sectors
"""

import shelve
import numpy as np

database = 'Outputs/BAEMC_set1' #file containig output from SCA analysis

#load SCA results
db_in = shelve.open(database)
D_seq = db_in['sequence']
D_sca = db_in['sca']
D_sector = db_in['sector']
db_in.close()

sectors = D_sector['ics']
Csca = D_sca['Csca']
ats = D_seq['ats']

avg_couplings = []
max_couplings = []
median_couplings =[]

#calculate coupling between sectors

i=0
for s1 in sectors:
    print i
    j = 0
    avg_couplings.append([])
    max_couplings.append([])
    median_couplings.append([])
    for s2 in sectors:
        if j <i:
            j+=1
            continue
        outfile_current_sector_couplings = open('outputs/coupling_between_sectors/between_sector_coupling_avg_'+ str(i+1) + '_' + str(j+1) +'.txt', 'w')
        outfile_current_sector_highly_coupled_pairs_selection = open('outputs/coupling_between_sectors/highly_coupled_positions_selection_'+ str(i+1) + '_' + str(j+1) +'.txt', 'w')
        outfile_current_sector_highly_coupled_pairs = open('outputs/coupling_between_sectors/highly_coupled_positions_'+ str(i+1) + '_' + str(j+1) +'.txt', 'w')
        highly_coupled_threshold= .5
        highly_coupled= []
        max_coupling = 0.0
        avg_coupling = 0.0
        total_pairs = 0.0
        all_couplings = []
        for p1 in s1.items:
            for p2 in s2.items:
                coupling = Csca[p1][p2]
                avg_coupling += coupling
                total_pairs += 1
                if coupling > max_coupling:
                    max_coupling = coupling
                if coupling > highly_coupled_threshold:
                    if p1 not in highly_coupled:
                        highly_coupled.append(p1)
                    if p2 not in highly_coupled:
                        highly_coupled.append(p2)
                    outfile_current_sector_highly_coupled_pairs.write(ats[p1] +", " + ats[p2] + "\n")
                outfile_current_sector_couplings.write(str(coupling) + "\n")
                all_couplings.append(coupling)
        for p in highly_coupled:
            outfile_current_sector_highly_coupled_pairs_selection.write("resi " + ats[p] + " + ")
        outfile_current_sector_highly_coupled_pairs_selection.close()
        outfile_current_sector_highly_coupled_pairs.close()
        outfile_current_sector_couplings.close()
        avg_couplings[i].append(avg_coupling/total_pairs)
        max_couplings[i].append(max_coupling)
        median_couplings[i].append(np.median(all_couplings))
        j +=1
    i+= 1
        

outfile_avg = open('outputs/coupling_between_sectors/between_sector_coupling_avg.txt', 'w')
outfile_max = open('outputs/coupling_between_sectors/between_sector_coupling_max.txt', 'w')
outfile_median = open('outputs/coupling_between_sectors/between_sector_coupling_median.txt', 'w')

#group sectors
groups  = {}
group_threshold =  .256 #mean + 2 standard deviations

for i in range(0, len(avg_couplings)):
    groups[i] = []
    for j in range(0, len(avg_couplings)):
        if j >= i:
            outfile_avg.write(str(avg_couplings[i][j-i]) + ",")
            outfile_max.write(str(max_couplings[i][j-i])+ ",")
            outfile_median.write(str(median_couplings[i][j-i])+ ",")
            if avg_couplings[i][j-i] >= group_threshold:
                groups[i].append(j+1)
        else:
            outfile_avg.write(str(avg_couplings[j][i-j])+ ",")
            outfile_max.write(str(max_couplings[j][i-j])+ ",")
            outfile_median.write(str(median_couplings[j][i-j])+ ",")
            if avg_couplings[j][i-j] >= group_threshold:
                groups[i].append(j+1)
    outfile_avg.write("\n")
    outfile_max.write("\n")
    outfile_median.write("\n")
outfile_avg.close()
outfile_max.close()
outfile_median.close()

for i in range(0, len(avg_couplings)):
    print "sector: " + str(i+1)
    print groups[i]
