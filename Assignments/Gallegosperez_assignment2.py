from Bio import SeqIO


#using the seq for the ref genome and store it as a string

for record in SeqIO.parse("Assignment2_refgenome.fasta", "fasta"):
    seq = record.seq


s_array = sorted([seq[i:] for i in range(len(seq))])

s_array[3]


    










