# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 12:03:45 2018

@author: Allison Walker
Script loads Hiseq data and outputs a list of mutations found in each sequencing read
"""

import os, os.path
import math
import operator
input_directory = "HiSeq/2_13_18" #directory containing Hiseq data
forward_primer_end = "GCAAGAC"
reverse_primer_end = "AAGGGTAT"
expected_reverse = "TCCTGGTCGGACATCAGGAGGTTAGTGCAATGGCATAAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTAT"
expected_forward = "GCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACTGAACATTGAGCCTTGATGTGTAGGATAGGTGGGAGGCTTTGAAGTGTGGACGCCAGTCTGCATGGAGCCGACCTTGAAATACCACCCTTTAATGTTTGATGTTCTAACGTTGACCCGTAATCCGGGTTGCGGACAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCTCCTCCTAAAGAGTAACGGAGGAGCACGAAGGTTGGCTAATCCTGGTCGGACATCAGGAGGTTAGTGCAATGGCATAAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTA"

quality_cutoff = 30
mutant_count_cutoff = 10
mutated_ranges = [(86, 98), (457, 481)]

"""
Takes the reverse complement of a DNA sequence

arguments:
    sequence - the sequence to take the reverse complement of
returns:
    the reverse complement of sequence
"""
def getReverseComplement(sequence):
    complement = ""
    for c in sequence:
        if c == 'A':
            complement += 'T'
        if c == 'T':
            complement += 'A'
        if c == 'G':
            complement += 'C'
        if c == 'C':
            complement += 'G'
        if c == 'N':
            complement += 'N'
    return getReverse(complement)

"""
reverses the order of a sequence

arguments:
    sequence - the sequence
returns:
    the sequence in reverse order
"""
def getReverse(sequence):
    reverse = ""
    for i in range(0, len(sequence)):
        reverse += sequence[len(sequence)-(i+1)]
    return reverse

"""
Compares a sequence to the expected sequence
arguments:
    sequence - the sequence to be compared to the expected sequence
    quality_sequence - a string containign the illumina quality scores for the sequence
    reverse - boolean, True if the sequence is a read in the reverse direction, False if it is a read in the forward direction

returns:
    (mutations, mutation_quality, mutation_to)
    mutations - a list of positions that are different from the expected sequence
    mutation_quality - quality at the mutated position
    mutation_to - a list of what the positions have been mutated to
"""
def compareToExpected(sequence, quality_sequence, reverse):
    mutations = []
    mutation_quality = []
    mutation_to = []
    if reverse:
        expected = expected_reverse
    else:
        expected = expected_forward
    for i in range(0, len(sequence)):
        if i > len(sequence)-1 or i > len(expected) -1:
            return (mutations, mutation_quality, mutation_to)
        if reverse != True and sequence[i] != expected[i]:
            if i > 230:
                continue
            mutations.append(i+2049)
            mutation_quality.append(quality_sequence[i])
            mutation_to.append(sequence[i])
        if reverse == True and sequence[len(sequence)-(i+1)] != expected[len(expected)-(i+1)]:
            if i > 240:
                continue
            mutations.append(2548-i)
            mutation_quality.append(quality_sequence[len(sequence)-(i+1)])
            mutation_to.append(sequence[len(sequence)-(i+1)])
    return (mutations, mutation_quality, mutation_to)

"""
Calculates the average quality at in a sequence read, also keeps track of the total quality score at each position
arguments:
    quality_at_each_pos - list that contains the total quality at each position
    quality_sequence - a sequence of the quality scores for a sequencing read
    num_sequences_with_position - a list containing the total number of sequences that contain this position after processing

"""
def calcAvgQuality(quality_at_each_pos, quality_sequence, num_sequences_with_position):
    avg=0.0
    if len(quality_at_each_pos)== 0:
        for q in quality_sequence:
            avg += ord(q)
            quality_at_each_pos.append(ord(q)-33)
            num_sequences_with_position.append(1.0)
    else:
        i=0
        for q in quality_sequence:
            avg += ord(q)
            if i < len(quality_at_each_pos):
                quality_at_each_pos[i]+=(ord(q)-33)
                num_sequences_with_position[i] += 1.0
            else:
                quality_at_each_pos.append(ord(q)-33)
                num_sequences_with_position.append(1.0)
            i+= 1
    avg /= len(quality_sequence)
    return (avg-33, quality_at_each_pos, num_sequences_with_position)

"""
Reads sequences from a file and writes a list of mutations in that sequence
arguments:
    fullpath - the path for the file containign the input sequence
    reverse - True if the sequence is a reverse read
    primer_end - primer used to amplify the DNA, used to recognize the start of the sequence of interest
"""
def getSequencesFromFile(fullpath, reverse, primer_end):
    input_file = open(fullpath, 'r')
    mutations_file = open(fullpath[0:fullpath.index(".")] + "_mutations.txt", 'w')
    avg_quality_file = open(fullpath[0:fullpath.index(".")] + "_avg_seq_quality.txt", 'w')
    avg_quality_file_kept_sequences = open(fullpath[0:fullpath.index(".")] + "_avg_quality_kept_sequences.txt", 'w')
    avg_quality_by_position = open(fullpath[0:fullpath.index(".")] + "_avg_quality_by_position.txt", 'w')
    quality_at_each_pos = []
    num_sequences_with_position= []
    line_index = 0
    total_sequences = 0
    total_sequences_with_primer = 0
    total_sequences_kept = 0
    index = ""
    sequence = ""
    quality_sequence = ""
    for line in input_file:
        if line_index%4 == 0: #a index line
            index = (line.replace("\n", ""))
        elif line_index%4 == 1: #a sequence line
            if reverse:
                sequence = (getReverseComplement(line.replace("\n", "")))
            else:
                sequence = (line.replace("\n", ""))
        elif line_index%4 == 2: #+ line
            line_index +=1
            continue
        elif line_index%4 == 3: #quality line
            quality_sequence = (line.replace("\n", ""))
            total_sequences += 1
            if primer_end in sequence:
                total_sequences_with_primer += 1
                if reverse:
                    processed_seq = sequence[0:sequence.index(primer_end)+len(primer_end)]
                    quality_sequence = quality_sequence[0:sequence.index(primer_end)+len(primer_end)]
                else:
                    processed_seq = sequence[sequence.index(primer_end):len(sequence)]
                    quality_sequence = quality_sequence[sequence.index(primer_end):len(sequence)]
                (mutations, mutation_quality, mutation_to) = compareToExpected(processed_seq, quality_sequence, reverse)
                (avg_quality, quality_at_each_pos_temp, num_sequences_with_position_temp) = calcAvgQuality(quality_at_each_pos, quality_sequence, num_sequences_with_position)
                avg_quality_file.write(str(avg_quality) +"\n")
                
                #skip reads with more than 30 mutations or an average quality score less than 20
                if len(mutations) > 30:
                    line_index +=1
                    continue
                if avg_quality <= 20:
                    line_index +=1
                    continue
                total_sequences_kept+=1
                (avg_quality, quality_at_each_pos, num_sequences_with_position) = calcAvgQuality(quality_at_each_pos, quality_sequence, num_sequences_with_position)
                avg_quality_file_kept_sequences.write(str(avg_quality) +"\n")
                mutations_file.write(index[0:43] + ":" + str(line_index) + "\n")
                for m in mutations:
                    mutations_file.write(str(m) + ",")
                mutations_file.write("\n")
                
                
                for q in mutation_quality:
                    if q != ",":
                        mutations_file.write(str(q) + ",")
                    else:
                        mutations_file.write("+,")
                mutations_file.write("\n")
                for m in mutation_to:
                    mutations_file.write(m + ",")
                mutations_file.write("\n")
        line_index +=1
        
    for i in range(0, len(quality_at_each_pos)):
        avg_quality_by_position.write(str(quality_at_each_pos[i]/num_sequences_with_position[i]) + "\n")
    sequence_count_file = open(fullpath[0:fullpath.index(".")] + "_sequence_counts.txt", 'w')
    sequence_count_file.write("total sequences: " + str(total_sequences) + "\n")
    sequence_count_file.write("total sequences with primer: " + str(total_sequences_with_primer) + "\n")
    sequence_count_file.write("total sequences kept: " + str(total_sequences_kept) + "\n")
    sequence_count_file.close()
    input_file.close()
    avg_quality_file.close()
    avg_quality_file_kept_sequences.close()
    avg_quality_by_position.close()
    mutations_file.close()

"""
Trims sequences to start at the beginning of the primer used to amplify the DNA
arguments:
    sequence_list - a list of sequences to be trimmed
    primer_end - primer sequence, the sequences are trimmed to start at this primer
    reverse - True if the read is a reverse read
    
returns:
    trimmed_list - a list of the trimmed sequences
"""
def trimSequences(sequence_list, primer_end, reverse):
    trimmed_list = []
    for sequence in sequence_list:
        if primer_end not in sequence or "N" in sequence:
            trimmed_list.append("")  #empty string if sequence doesn't have primer sequence
            continue
        else:
            if reverse:
                trimmed_list.append(sequence[0:sequence.index(primer_end)+len(primer_end)])
            else:
                trimmed_list.append(sequence[sequence.index(primer_end):len(sequence)])
    return trimmed_list

"""
Joins forward and reverse reads
"""
def joinEnds(r_seq, f_seq, q_seq_f, q_seq_r, expected_sequence):
    exact_match_f = ""
    exact_match_r = ""
    exact_match_q_f = ""
    exact_match_q_r = ""
    for i in range(0, len(f_seq)):
        if f_seq[i] == expected_sequence[i]:
            exact_match_f += f_seq[i]
            exact_match_q_f += q_seq_f[i]
        elif expected_sequence[i] == "N":
            exact_match_f += f_seq[i]
            exact_match_q_f += q_seq_f[i]
        else:
            break
    for i in range(0, len(r_seq)):
        if r_seq[len(r_seq)-i - 1] == expected_sequence[len(expected_sequence)-i - 1]:
            exact_match_r += r_seq[len(r_seq)-i - 1]
            exact_match_q_r += q_seq_r[len(r_seq)-i - 1]
        elif expected_sequence[len(expected_sequence)-i - 1] == "N":
            exact_match_r += r_seq[len(r_seq)-i - 1]
            exact_match_q_r += q_seq_r[len(r_seq)-i - 1]
        else:
            break

    exact_match_r = getReverse(exact_match_r)
    exact_match_q_r = getReverse(exact_match_q_r)
    overlap = len(exact_match_r) + len(exact_match_f) - len(expected_sequence) #length of sequence

    if overlap  >= 0:
        return (exact_match_f[0:len(exact_match_f)-overlap] + exact_match_r, exact_match_q_f[0:len(exact_match_f)-overlap] + exact_match_q_r)
    else:
        return ("", "")
        if overlap < -80 or len(exact_match_f) < 30 or len(exact_match_r) < 50:
            return ""
        elif "N" in exact_match_f:
            return ""
        elif "N" in exact_match_r:
            return ""
        else:
            gap_Ns = ""
            for i in range(0, -1*overlap):
                gap_Ns += "N"

            return exact_match_f + gap_Ns + exact_match_r

"""
Determines if a position is in the parts of the sequence where mutations were made
"""        
def inMutatedRange(index, mutated_ranges):
    for r in mutated_ranges:
        if index >= r[0] or index <= r[1]:
            return True
    return False

"""
Joins forward and reverse reads
"""
def joinEndsNoExactMatch(r_seq, f_seq, q_seq_f, q_seq_r, expected_sequence):
    mismatch_positions =[]
    mismatch_positions =[]
    for i in range(0, len(f_seq)):
        if f_seq[i] != expected_sequence[i]:
            #check quality of mismatch
            mismatch_q = q_seq_f[i]
            if ord(mismatch_q) -33 > quality_cutoff:# or inMutatedRange(i, mutated_ranges):
                mismatch_positions.append(i)
    for i in range(0, len(r_seq)):
        if r_seq[len(r_seq)-i - 1] != expected_sequence[len(expected_sequence)-i - 1]:
            mismatch_q = q_seq_r[len(r_seq)-i - 1]
            if ord(mismatch_q) -33 > quality_cutoff or inMutatedRange(len(expected_sequence) -i-1, mutated_ranges):
                mismatch_positions.append(len(expected_sequence) -i-1)

    overlap = len(r_seq) + len(f_seq) - len(expected_sequence) #length of sequence
    if overlap  >= 0:
        return (f_seq[0:len(exact_match_f)-overlap] + r_seq, q_seq_f[0:len(exact_match_f)-overlap] + q_seq_r, mismatch_positions)
    num_Ns = abs(overlap)
    N_seq = ""
    for i in range(0, num_Ns):
        N_seq += "N"
    return (f_seq + N_seq + r_seq, q_seq_f + N_seq + q_seq_r, mismatch_positions)


"""
Gets the counts for each mutant
"""
def getMutCounts(sequences, quality_seqs, mutant_positions):
    counts = {}
    i = 0
    for s in sequences:
        q = quality_seqs[i]
        if s == "" :
            i += 1
            continue


        mutant_seq = ""
        qualities = []
        quality_too_low = False
        for p in mutant_positions:
            mutant_seq += s[p]
            qualities.append((ord(q[p]) - 33))
            if (ord(q[p]) - 33) < quality_cutoff:
               quality_too_low = True
        if quality_too_low:
            i +=1
            continue
        if mutant_seq not in counts:
            counts[mutant_seq] = 1
        else:
            counts[mutant_seq] += 1
        i +=1
    return counts


forward_sequences = {}
reverse_sequences = {}
quality_seqs_f = {}
quality_seqs_r = {}
sample_list = []
reverse = False
mutations_f = {}
mutations_r = {}
sample_names = []
#read in sequences from file and trim sequences based on primer seqs
print input_directory
for root, _, files in os.walk(input_directory):
    for f in files:
        fullpath = os.path.join(root, f)
        print fullpath
        if "zip" in fullpath or "gz" in fullpath or "mutations" in fullpath or "_processed" in fullpath or "mutant" in fullpath or "html" in fullpath or "fastqc" in fullpath or "quality" in fullpath:
            continue
        if "data" in fullpath:
            continue
        if "sequence_counts" in fullpath:
            continue
        
        #if "E1S8" not in fullpath:
         #   continue
        #if "E1S1" in fullpath or "E1S2" in fullpath or "E1S3" in fullpath:
         #   continue
        #print fullpath
        if 'R1' in fullpath:
            reverse = True
            primer_end = reverse_primer_end
        else:
            reverse = False
            primer_end = forward_primer_end
        if "_" not in f:
            continue
        sample_name = f[0:f.index("_")]
        print sample_name
        sample_list.append(sample_name)
        if sample_name not in sample_names:
            sample_names.append(sample_name)
        mutations = getSequencesFromFile(fullpath, reverse, primer_end)
        if reverse:
            mutations_r[sample_name] = mutations
        else:
            mutations_f[sample_name] = mutations

#remove low quality or unknown mutations        
min_count_threshold = 2
orig_mutant_counts = {}
total_counts={}

sequence_counts = {}
for sn in sample_names:
    sequence_counts[sn] = 0      
    
for sample in sample_names:
    if "E1" in sample:
        lane_numer =1
    else:
        lane_number = 2
    forward_file = open(input_directory + "/" + sample + "_L00"+str(lane_number)+"_R2_001_mutations.txt", 'r')
    reverse_file = open(input_directory + "/" + sample + "_L00"+str(lane_number)+"_R1_001_mutations.txt", 'r')
    WT_count = 0
    print(sample)

    mutant_counts = {}
    while True:
        f_line = forward_file.readline()
        if not f_line:
            break
        r_line = reverse_file.readline()
        if not r_line:
            break
        while f_line != r_line and f_line and r_line:
            if len(f_line.split(":")) < 7:
                continue
            if len(r_line.split(":")) < 7:
                continue
            if int(f_line.split(":")[7]) < int(r_line.split(":")[7]):
                while f_line and int(f_line.split(":")[7]) < int(r_line.split(":")[7]):
                    for i in range(0, 4):
                        f_line = forward_file.readline()
                    if len(f_line.split(":")) < 7:
                        break
            else:
                while r_line and int(f_line.split(":")[7]) > int(r_line.split(":")[7]) and r_line:
                    for i in range(0, 4):
                        r_line = reverse_file.readline()
                    if len(r_line.split(":")) < 7:
                        break
        if len(f_line.split(":"))<7 or len(r_line.split(":"))< 7:
            continue
        if int(f_line.split(":")[7]) != int(r_line.split(":")[7]):
            print "LINE INDEX DOESN'T MATCH!"
            break
        forward_mut = forward_file.readline().replace("\n", "").split(",")
        forward_mut = forward_mut[0:len(forward_mut)-1]
        forward_quality = forward_file.readline().replace("\n", "").split(",")
        forward_quality = forward_quality[0:len(forward_quality)-1]
        forward_mut_to = forward_file.readline().replace("\n", "").split(",")
        reverse_mut = reverse_file.readline().replace("\n", "").split(",")
        reverse_mut = reverse_mut[0:len(reverse_mut)-1]
        reverse_quality = reverse_file.readline().replace("\n", "").split(",")
        reverse_quality = reverse_quality[0:len(reverse_quality)-1]
        reverse_mut_to = reverse_file.readline().replace("\n", "").split(",")
        
        
        sequence_counts[sample] += 1
        mut_key = ""
        if len(forward_mut) > len(forward_mut_to):
            continue
        for i in range(0, len(forward_mut)):
            if int(forward_mut[i][0:len(forward_mut[i])]) > 2256:
                continue
            if ord(forward_quality[i]) - 33 >20 or orig_mutant_counts[sample][forward_mut[i] + forward_mut_to[i]]>min_count_threshold:#.00001*total_counts[sample]:
                mut_key += forward_mut[i] + forward_mut_to[i] + " "

        for i in range(0, len(reverse_mut)):
            if int(reverse_mut[i][0:len(reverse_mut[i])]) < 2400:
                continue
            if ord(reverse_quality[i]) - 33 >20 or orig_mutant_counts[sample][reverse_mut[i] + reverse_mut_to[i]]>min_count_threshold:#.00001*total_counts[sample]:
                mut_key += reverse_mut[i] + reverse_mut_to[i] + " "
        if mut_key not in mutant_counts:
            mutant_counts[mut_key] = 1
        else:
            mutant_counts[mut_key] += 1
        if mut_key == "":
            WT_count += 1
        
    outfile = open(input_directory + "/" + sample + "_mutant_counts_quality_threshold.txt", 'w')
    outfile.write("WT count, " + str(WT_count) + "\n")
    for mut_key in mutant_counts:
        if mutant_counts[mut_key] < 5:
            continue
        outfile.write(mut_key + ", " + str(mutant_counts[mut_key]) + "\n")
    outfile.close()
print(sequence_counts)