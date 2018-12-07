# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 18:00:08 2018

@author: allis
"""

import math
import numpy as np
import shelve
import matplotlib as plt
import scipy
from scipy.stats.stats import pearsonr   

groups = [[0, 2, 4, 12], [1,22],[3,13],[0, 1, 2, 3, 4,12,13,22]] #grouped ICs into sectors 1, 2, 3, and all sectors

database = 'Outputs/BAEMC_set1'

outputfile_root = "outputs/fixed_sector_position_ranking_sector_"

#load SCA results
db_in = shelve.open(database)
D_seq = db_in['sequence']
D_sca = db_in['sca']
D_sector = db_in['sector']
db_in.close()
ats = D_seq['ats']
sectors = D_sector['ics']
Vp_eig = D_sector['Vsca']
Wpica = D_sector['Wpica']
Vp = D_sector['Vpica']
eigen_vals = D_sector['Lsca']
Csca = D_sca['Csca']
icsize =D_sector['icsize']
sortedpos = D_sector['sortedpos']

#calculate weights for each IC
ev_weights = []
for i in range(0, Vp.shape[1]):
    ev= np.matmul(Vp[:,i].T, np.matmul(Csca, Vp[:,i]))
    ev_weights.append(ev)
    
i= 0
for group in groups:
    i += 1
    outfile = open(outputfile_root   + str(i) + ".txt", 'w')
    average_vector = np.zeros(Vp[:,0].shape)
    total_weights = 0
    pos_in_sector = []
    for j in group:
        average_vector += ev_weights[j]*Vp[:,j]
        total_weights += ev_weights[j]
        current_size = icsize[j]
        previous_sizes = int(np.sum(icsize[0:j]))
        pos_in_sector += sortedpos[previous_sizes:previous_sizes + current_size]
    average_vector /= total_weights
    average_vector = average_vector[sorted(pos_in_sector)]
    sorted_indices = np.argsort(-1*average_vector)
    for index in sorted_indices:
        outfile.write(ats[sorted(pos_in_sector)[index]] + ", " + str(average_vector[index]) +  "\n")
    outfile.close()