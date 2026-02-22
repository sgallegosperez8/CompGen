import bisect
from Bio import SeqIO


#using the seq for the ref genome and store it as a string

for record in SeqIO.parse("Assignment2_refgenome.fasta", "fasta"):
    seq = record.seq


def suffixArray(s):
    suffixes = [(s[i:], i) for i in range(len(s))]
    suffixes.sort(key=lambda x: x[0])
    suffixes
    return [s[1] for s in suffixes]

suffixa=suffixArray('banana')
len(suffixa)
sorted(suffixa)


target = 'na'
i=0

#creating while loop until false
while i < len(suffixa):
    print(suffixa[i])
    i +=1





    










