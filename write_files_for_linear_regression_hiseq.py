# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 15:19:24 2018

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
experiment_index = 2
input_directory = "HiSeq/2_13_18" #directory where processed sequence reads are stored
output_directory = "HiSeq/2_13_18/linear_regression_data/" #directory to write linear regression to
thresholded_slopes_out = open("HiSeq/2_13_18/final_data/HiSeq_mutant_matrix_slope_18_26_"+str(experiment_index)+"_"+ str(run_index)+"_threshold" + str(mutant_count_treshold) + "_" + "r2_" + str(r2_threshold) +".txt", 'w')

sample_indices = {}
#E1 run 1
if experiment_index == 1 and run_index == 2:
    sample_indices[(2, 0)] ="E1S1_ATTACTCG-TATAGCCT"
    sample_indices[(3, 0)] ="E1S2_ATTACTCG-ATAGAGGC"
    sample_indices[(run_index, 18)] ="E1S3_ATTACTCG-GGCTCTGA"
    sample_indices[(run_index, 20)] ="E1S5_ATTACTCG-GTACTGAC"
    sample_indices[(run_index, 24)] ="E1S7_TCCGGAGA-ATAGAGGC"
    sample_indices[(run_index, 26)] ="E1S9_TCCGGAGA-GGCTCTGA"
#E1 run 2
elif experiment_index == 1 and run_index == 3:
    sample_indices[(2, 0)] ="E1S1_ATTACTCG-TATAGCCT"
    sample_indices[(3, 0)] ="E1S2_ATTACTCG-ATAGAGGC"
    sample_indices[(run_index, 18)] ="E1S4_ATTACTCG-CAGGACGT"
    sample_indices[(run_index, 20)] ="E1S6_TCCGGAGA-TATAGCCT"
    sample_indices[(run_index, 24)] ="E1S8_TCCGGAGA-CCTATCCT"
    sample_indices[(run_index, 26)] ="E1S10_TCCGGAGA-CAGGACGT"
#E2 run 1
elif experiment_index == 2 and run_index == 2:
    sample_indices[(2, 0)] ="E2S1_ATTACTCG-TATAGCCT"
    sample_indices[(3, 0)] ="E2S2_ATTACTCG-ATAGAGGC"
    sample_indices[(run_index, 18)] ="E2S3_ATTACTCG-CCTATCCT"
    sample_indices[(run_index, 20)] ="E2S5_ATTACTCG-CAGGACGT"
    sample_indices[(run_index, 24)] ="E2S7_TCCGGAGA-TATAGCCT"
    sample_indices[(run_index, 26)] ="E2S9_TCCGGAGA-CCTATCCT"
#E2 run 2
elif experiment_index== 2 and run_index==3:
    sample_indices[(2, 0)] ="E2S1_ATTACTCG-TATAGCCT"
    sample_indices[(3, 0)] ="E2S2_ATTACTCG-ATAGAGGC"
    sample_indices[(run_index, 18)] ="E2S4_ATTACTCG-GGCTCTGA"
    sample_indices[(run_index, 20)] ="E2S6_ATTACTCG-GTACTGAC"
    sample_indices[(run_index, 24)] ="E2S8_TCCGGAGA-ATAGAGGC"
    sample_indices[(run_index, 26)] ="E2S10_TCCGGAGA-GGCTCTGA"
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
            for m in mutant.split(" "):
                if int(m[0:len(m)-1]) > 2300 and int(m[0:len(m)-1]) < 2400:
                    continue
            if (int(line.split(",")[1].replace("\n", ""))) > mutant_count_treshold:
                if mutant not in mutants:
                    mutants[mutant] = {}
                mutants[mutant][k] = int(line.split(",")[1].replace("\n", ""))
                
print mutants             
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
    if experiment_index == 1:
        time_points = [17.5,20,24,26]
    else:
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
