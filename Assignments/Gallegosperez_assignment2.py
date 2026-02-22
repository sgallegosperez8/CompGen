import bisect
import re
from Bio import SeqIO

#open the read to be searched in the ref genome
file = open("Assignment2_read.fasta","r")

#using the seq for the ref genome and store it as a string
for record in SeqIO.parse("Assignment2_refgenome.fasta", "fasta"):
    seq = record.seq

#run it through function suffixArray() to create suffix array/tree of reference genome 
def suffixArray(s):
    suffix = [(str(s[i:]), i) for i in range(len(s))]
    suffix.sort(key=lambda x: x[0])
    return suffix

suffixa=suffixArray(seq)

#readseq() to get the target sequence
def readseq(file):
    pat=r"[ATCG]"
    for line in file:
        if re.search(pat,line):
            tar=line.strip()
    return tar

tar = readseq(file)

#b_search() to use binary search to locate all positions 
def b_search(tar, suffix):
    mini = 0
    max_s = (len(suffix))-1
    match=[]
    while mini<=max_s:
        mid = mini + (max_s-mini)//2
        c_str = suffix[mid][0]
        if re.match(tar,suffix[mid][0]):
            match.append(suffix[mid][1])
            for_val=mid+1
            rev_val=mid-1
            while rev_val >=0 and re.match(tar,suffix[rev_val][0]):
                match.append(suffix[rev_val][1])
                rev_val -=1

            while rev_val < len(suffix) and re.match(tar,suffix[for_val][0]):
                match.append(suffix[for_val][1])
                for_val +=1
            break

        if c_str < tar:
            mini = mid +1
        else:
            max_s= mid -1
    return match


match=b_search(tar,suffixa) 


print(f"Positions for {tar}: {match}")
