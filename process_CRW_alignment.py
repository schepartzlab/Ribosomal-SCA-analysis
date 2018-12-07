# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 18:04:07 2018

@author: Allison Walker

script process CRW 23S alignments for bacteria, archaea, eukaryota, mitochondria, and chloroplasts
CRW alignment files should be saved in the same folder as this script with the filenames:
    23S.A.alnfasta for achaeal sequences
    23S.B.ALL.alnfasta for bacterial sequences
    23S.E.alnfasta for eukaryotic sequences
    23S.C.alnfasta for chloroplast sequences
    23S.M.alnfasta for mitochondrial sequences
writes two new alignment files:
    CRW_seqs_BAEMC.fasta which contains all sequences, with gaps in the E. coli sequence
    CRW_no_gaps_BAEMC.fasta contains all sequences but with all positions that are gaps in the E. coli reference sequence removed, this alignment was used
    as the input for SCA analysis
"""

"""
Functionn takes a CRW alignment file, identifies the E. coli reference sequence and returns
a list of sequences in the file.

Arguments:
    input_file the file containig the CRW alignmetn
    start_seq expected sequence at the start of the e coli reference sequence
    end_seq expected sequence at the end of the e coli reference squence
    
Returns: (new_titles, sequences, e_coli_seq[start_index:end_index])
    new_titles a list of the headers for each sequence
    sequences a list of sequences
    e_coli_seq[start_index:end_index] the e coli reference sequence
    
"""
def processFile(input_file, start_seq, end_seq):
    titles = []
    sequences = []
    seq = ""
    e_coli_seq = ""
    title = ""
    for line in input_file:
        if ">" in line:
            if "Escherichia_coli" in title:
                e_coli_seq = seq
            if "A" in seq.upper() and len(seq) > 0 and len(title) > 0:
                titles.append(title)
                sequences.append(seq.upper())
                seq = ""           
            title = line

        else:
            seq += line.replace("\n", "")
    titles.append(title)
    sequences.append(seq)
    
    start_index = e_coli_seq.find(start_seq)
    end_index = e_coli_seq.find(end_seq) + len(end_seq)
    new_sequences = []
    for s in sequences:
        new_sequences.append(s[start_index:end_index])

    to_add = []
    new_titles = []
    for i in range(0, len(new_sequences)):
        if 1.0*new_sequences[i].count("-")/len(new_sequences[i]) < .7:
            to_add.append(i);

    
    sequences = [];
    for i in range(0, len(new_sequences)):
        if i in to_add:
            sequences.append(new_sequences[i])
            new_titles.append(titles[i])

    return (new_titles, sequences, e_coli_seq[start_index:end_index])

#open output files
out = open("CRW_seqs_BAEMC.fasta", "w")
out_no_gaps = open("CRW_no_gaps_BAEMC.fasta", "w")

#open and process all CRW files
arch = open("23S.A.alnfasta", 'r')
(arch_titles, arch_seqs, a_coli) = processFile(arch, "GGUU-AAGC-", "GCUUAA-CC-UU")
arch.close()

bac = open("23S.B.ALL.alnfasta", 'r')
(bac_titles, bac_seqs, b_coli) = processFile(bac, "GGUU-AAGC-", "GCUUAA-CC-UU")
bac.close()

chlor = open("23S.C.alnfasta", 'r')
(chlor_titles, chlor_seqs, c_coli) = processFile(chlor, "GGUUAAGCG-AC", "GCUUAACC-UU")
chlor.close()

mito = open("23S.M.alnfasta", 'r')
(mito_titles, mito_seqs, m_coli) = processFile(mito, "GGU-UAAGC", "GCUUAACCUU")
mito.close()

euk = open("23S.E.alnfasta", 'r')
(euk_titles, euk_seqs, e_coli) = processFile(euk, "GGUUAAGC", "GCUUAACCUU");
euk.close()

#line up the gaps
new_a = [];
new_a_ng = [];
for i in range(0, len(arch_seqs)):
    new_a.append("");
    new_a_ng.append("");
new_b = [];
new_b_ng = [];
for i in range(0, len(bac_seqs)):
    new_b.append("");
    new_b_ng.append("");
new_c = [];
new_c_ng = [];
for i in range(0, len(chlor_seqs)):
    new_c.append("");
    new_c_ng.append("")
new_m = [];
new_m_ng = [];
for i in range(0, len(mito_seqs)):
    new_m.append("")
    new_m_ng.append("")
new_e = [];
new_e_ng = [];
for i in range(0, len(euk_seqs)):
    new_e.append("")
    new_e_ng.append("")
i_a = 0;
i_b = 0;
i_c = 0;
i_m = 0;
i_e = 0;
new_coli = "";
while i_a < len(a_coli) and i_b < len(b_coli) and i_c < len(c_coli) and i_m < len(m_coli) and i_e < len(e_coli):
    a_c = (a_coli[i_a] == "-");
    b_c = (b_coli[i_b] == "-");
    c_c = (c_coli[i_c] == "-");
    m_c = (m_coli[i_m] == "-");
    e_c = (e_coli[i_e] == "-");
    if not (a_c or b_c or c_c or m_c or e_c):
        #don't need to add any dashes
        for j in range(0, len(arch_seqs)):
            new_a[j] = new_a[j] + arch_seqs[j][i_a]
        for j in range(0, len(bac_seqs)):
            new_b[j] = new_b[j] + bac_seqs[j][i_b]
        for j in range(0, len(chlor_seqs)):
            new_c[j] = new_c[j] + chlor_seqs[j][i_c]
        for j in range(0, len(mito_seqs)):
            new_m[j] = new_m[j] + mito_seqs[j][i_m]
        for j in range(0, len(euk_seqs)):
            new_e[j] = new_e[j] + euk_seqs[j][i_e]
        new_coli = new_coli + a_coli[i_a];
        i_a +=1;
        i_b +=1;
        i_c +=1;
        i_m +=1;
        i_e +=1;
    else:
        if a_c:
             for j in range(0, len(arch_seqs)):
                new_a[j] = new_a[j] + arch_seqs[j][i_a]
             i_a +=1;
             new_coli = new_coli + a_coli[i_a];
        else:
            for j in range(0, len(arch_seqs)):
                new_a[j] = new_a[j] +"-"
            new_coli = new_coli + "-"

        if b_c:
             for j in range(0, len(bac_seqs)):
                new_b[j] = new_b[j] + bac_seqs[j][i_b]
             i_b +=1;
        else:
            for j in range(0, len(bac_seqs)):
                new_b[j] = new_b[j] +"-"

        if c_c:
             for j in range(0, len(chlor_seqs)):
                new_c[j] = new_c[j] + chlor_seqs[j][i_c]
             i_c +=1;
        else:
            for j in range(0, len(chlor_seqs)):
                new_c[j] = new_c[j] +"-"

        if m_c:
             for j in range(0, len(mito_seqs)):
                new_m[j] = new_m[j] + mito_seqs[j][i_m]
             i_m +=1;
        else:
            for j in range(0, len(mito_seqs)):
                new_m[j] = new_m[j] +"-"

        if e_c:
             for j in range(0, len(euk_seqs)):
                new_e[j] = new_e[j] + euk_seqs[j][i_e]
             i_e +=1;
        else:
            for j in range(0, len(euk_seqs)):
                new_e[j] = new_e[j] +"-"

#write resulting alignments
for i in range(0, len(arch_seqs)):
    out.write(arch_titles[i]);
    out.write(new_a[i] + "\n");
for i in range(0, len(bac_seqs)):
    out.write(bac_titles[i]);
    out.write(new_b[i] + "\n");
for i in range(0, len(chlor_seqs)):
    out.write(chlor_titles[i]);
    out.write(new_c[i] + "\n");
for i in range(0, len(mito_seqs)):
    out.write(mito_titles[i]);
    out.write(new_m[i] + "\n");
for i in range(0, len(euk_seqs)):
    out.write(euk_titles[i]);
    out.write(new_e[i] + "\n");

new_coli = new_a[0];
for i in range(0, len(new_coli)):
    if new_coli[i] == "-":
        continue;
    else:
        for j in range(0, len(new_a)):
            new_a_ng[j] = new_a_ng[j] + new_a[j][i];
        for j in range(0, len(new_b)):
            new_b_ng[j] = new_b_ng[j] + new_b[j][i];
        for j in range(0, len(new_e)):
            new_e_ng[j] = new_e_ng[j] + new_e[j][i];
           # print len(new_e[j])
            #new_e_ng[j] = new_e_ng[j];
        for j in range(0, len(new_c)):
            new_c_ng[j] = new_c_ng[j] + new_c[j][i];
        for j in range(0, len(new_m)):
            new_m_ng[j] = new_m_ng[j] + new_m[j][i];
arch_count = 0
bac_count = 0
euk_count = 0
mito_count = 0
chlor_count = 0
for i in range(0, len(arch_seqs)):
    if "N" in arch_seqs[i]:
        continue
    if "coli" not in arch_titles[i]:
        arch_count += 1
    out_no_gaps.write(arch_titles[i]);
    out_no_gaps.write(new_a_ng[i] + "\n");
print "There are: " +str(arch_count) + " arch seqs"
for i in range(0, len(bac_seqs)):
    if "N" in bac_seqs[i]:
        continue
    if "REFERENCE" not in bac_titles[i]:
        bac_count += 1
    out_no_gaps.write(bac_titles[i]);
    out_no_gaps.write(new_b_ng[i] + "\n");
print "There are: " +str(bac_count) + " bac seqs"
for i in range(0, len(chlor_seqs)):
    if "N" in chlor_seqs[i]:
        continue
    if "REFERENCE" not in chlor_titles[i]:
        chlor_count += 1
    out_no_gaps.write(chlor_titles[i]);
    out_no_gaps.write(new_c_ng[i] + "\n");
print "There are: " +str(chlor_count) + " chlor seqs"
for i in range(0, len(mito_seqs)):
    if "N" in mito_seqs[i]:
        continue
    if "REFERENCE" not in mito_titles[i]:
        mito_count += 1
    out_no_gaps.write(mito_titles[i]);
    out_no_gaps.write(new_m_ng[i] + "\n");
print "There are: " +str(mito_count) + " mito seqs"
for i in range(0, len(euk_seqs)):
    if "N" in euk_seqs[i]:
        continue
    if "REFERENCE" not in euk_titles[i]:
        euk_count += 1
    out_no_gaps.write(euk_titles[i]);
    out_no_gaps.write(new_e_ng[i] + "\n");
print "There are: " +str(euk_count) + " euk seqs"
out.close();
out_no_gaps.close();

