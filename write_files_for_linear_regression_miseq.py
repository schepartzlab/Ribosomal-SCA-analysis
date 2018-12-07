# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 15:27:02 2018

@author: Allison Walker
Script peforms linear regression on time vs log(relative frequency) for Miseq experiments
"""

from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
from sklearn import datasets, linear_model
import math

nucleotides = ["A", "C", "G", "T"]

mutant_count_treshold = 300 #the required number of reads containing each mutant for the mutant to be included in the dataset
r2_threshold =.6 #the required R^2 on the growth curve in order for the mutant to be included in the dataset 
run_index=2
input_directory = "12_21_17_sequencing/171222-0437M_Link_Yang_AW1Miseq/fastq_Lane1"
output_directory = "12_21_17_sequencing/171222-0437M_Link_Yang_AW1Miseq/fastq_Lane1/linear_regression_data/threshold_" + str(mutant_count_treshold) + "/"
thresholded_slopes_out = open("12_21_17_sequencing/171222-0437M_Link_Yang_AW1Miseq/fastq_Lane1/final_data/12_21_17_mutant_matrix_slope_18_26_"+str(run_index)+"_threshold" + str(mutant_count_treshold) + "_" + "r2_" + str(r2_threshold) +".txt", 'w')
sample_indices = {}
sample_indices[(2, 0)] ="0-2_S1"
sample_indices[(3, 0)] ="0-3_S2"
if run_index == 2:
    sample_indices[(run_index, 18)] ="1-"+str(run_index)+"_S3"
    sample_indices[(run_index, 20)] ="2-"+str(run_index)+"_S5"
    sample_indices[(run_index, 24)] ="3-"+str(run_index)+"_S7"
    sample_indices[(run_index, 26)] ="4-"+str(run_index)+"_S9"
elif run_index == 3:
    sample_indices[(run_index, 18)] ="1-"+str(run_index)+"_S4"
    sample_indices[(run_index, 20)] ="2-"+str(run_index)+"_S6"
    sample_indices[(run_index, 24)] ="3-"+str(run_index)+"_S8"
    sample_indices[(run_index, 26)] ="4-"+str(run_index)+"_S10"


WT_seq = ""
nucleotides = ["A", "C", "G", "T"]
time_point_list= [(run_index, 18), (run_index, 20), (run_index, 24), (run_index, 26)]
mutant_matrix_slopes = {}
                  
mutants = {}
for k in sample_indices:
    sample_name = sample_indices[k]
    in_file = open(input_directory + "/" + sample_name + "_mutant_counts_quality_threshold.txt")
    for line in in_file:
        if "WT" in line:
            if "WT" not in mutants:
                mutants["WT"] = {}
            mutants["WT"][k] = int(line.split(",")[1].replace("\n", ""))
        else:
            if line.split(",")[0] == "":
                continue
            mutant = line.split(",")[0][0: len(line.split(",")[0])-1]
            if (int(line.split(",")[1].replace("\n", ""))) > mutant_count_treshold:
                if mutant not in mutants:
                    mutants[mutant] = {}
                mutants[mutant][k] = int(line.split(",")[1].replace("\n", ""))
                
               
for m in mutants:
    if (2,0) not in mutants[m] or (3,0) not in mutants[m]:
        continue
    if (run_index, 18) not in mutants[m]:
        continue
    if (run_index, 20) not in mutants[m]:
        continue
    if (run_index, 24) not in mutants[m]:
        continue
    if (run_index, 26) not in mutants[m]:
        continue
    if "WT" in m:
        continue
    avg_tp_zero_mut = (mutants[m][(2, 0)] + mutants[m][(3, 0)])/2.0
    avg_tp_zero_WT = (mutants["WT"][(2, 0)] + mutants["WT"][(3, 0)])/2.0
    relative_freq_18 = math.log(1.0*mutants[m][(run_index, 18)]/ avg_tp_zero_mut) - math.log(1.0*mutants["WT"][(run_index, 18)]/avg_tp_zero_WT)
    relative_freq_20 = math.log(1.0*mutants[m][(run_index, 20)]/ avg_tp_zero_mut) - math.log(1.0*mutants["WT"][(run_index, 20)]/avg_tp_zero_WT)
    relative_freq_24 = math.log(1.0*mutants[m][(run_index, 24)]/ avg_tp_zero_mut) - math.log(1.0*mutants["WT"][(run_index, 24)]/avg_tp_zero_WT)
    relative_freq_26 = math.log(1.0*mutants[m][(run_index, 26)]/ avg_tp_zero_mut) - math.log(1.0*mutants["WT"][(run_index, 26)]/avg_tp_zero_WT)
    enrichments_18_26 = [relative_freq_18,relative_freq_20, relative_freq_24, relative_freq_26]
    time_points = [18,20,24,26]
    regr = linear_model.LinearRegression()
    regr.fit(np.array(time_points).reshape(4,1), np.array(enrichments_18_26).reshape(4,1))
    enrichments_pred = regr.predict(np.array(time_points).reshape(4,1))
    growth_rate_1 = regr.coef_[0][0]
    r2_1 = r2_score(enrichments_18_26, enrichments_pred)
    if len(m.split(" ")) < 2 and "WT" not in m or len(m.split(" ")) == 2:
        outfile = open(output_directory + m + "_" + str(run_index) +".txt", 'w')
        outfile.write(str(growth_rate_1) + "," + str(r2_1) + "\n")
        outfile.write("18," + str(relative_freq_18) + "," + "\n")
        outfile.write("20," + str(relative_freq_20) + "," + "\n")
        outfile.write("24," + str(relative_freq_24) + "," + "\n")
        outfile.write("26," + str(relative_freq_26)+ ","  + "\n")
        outfile.close()
    if len(m.split(" ")) < 2:
        if r2_1 > r2_threshold/100.0:
            mut_pos = m.split(" ")[0][0:len(m.split(" ")[0])-1]
            mut_to = m.split(" ")[0][len(m.split(" ")[0])-1:len(m.split(" ")[0])]
            if mut_pos not in mutant_matrix_slopes:
                mutant_matrix_slopes[mut_pos] = {}
            mutant_matrix_slopes[mut_pos][mut_to] = growth_rate_1
print mutant_matrix_slopes
for loc in mutant_matrix_slopes:
   thresholded_slopes_out.write(str(loc) + ",") 
   for n in nucleotides:
       if n not in mutant_matrix_slopes[loc]:
            thresholded_slopes_out.write("0,")
       else:
            thresholded_slopes_out.write(str(mutant_matrix_slopes[loc][n]) +",")
   thresholded_slopes_out.write("\n")
thresholded_slopes_out.close()
        
    
