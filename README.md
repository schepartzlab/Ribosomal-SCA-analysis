# Ribosomal-SCA-analysis

Scripts for analyzing the output of SCA 6 on ribosomal RNA and for analyzing next generation sequencing data from continuous culture experiments

process_CRW_alignment.py combines alignments from the CRW database into one master alignment that can be used as the input for SCA analysis

group_sectors.py groups ICs from SCA analysis into sectors based on the coupling between residues in each ICs

rank_positions_by_weighted_contribution_avg.py calculates the weighted average of each posisitions contribution to each sector as well as all three 
sectors combined.

Hiseq_data_analysis.py and Miseq_data_analysis.py executes the intial steps of analyzing next generation sequencing data. They determine the mutations
present in each sequencing read

write_files_for_linear_regression_hiseq.py and write_files_for_linear_regression_miseq.py perform the final steps of analyzing next generaion sequencing
data. They use the output from Hiseq_data_analysis.py and Miseq_data_analysis.py. These scripts perform linear regressions on the log of the relative frequency
of each mutant over time.

environment.yml file contains information about all package versions used to produce the results described in the paper, although we anticipate that the same
results would be obtained using any installation of the Anaconda build of Python 2.

SCA_matrix.csv contains the SCA coupling matrix, all values have been rounded to three digits to reduce the file size. Therefore, this matrix should only be used for
visualization or ranking pairs of positions by degree of coupling, not to reproduce the ICs. To reproduce the ICs described in the paper please run the Jupyter notebook.

Jupyter notebook folder contains a jupyter notebook with a demonstration of the code used to produce the results described in the paper.