# Computational Genomic Projects 

## 1: [HMM BaseCaller and FastQC](https://github.com/sgallegosperez8/CompGen/tree/main/HMM-BaseCaller-and-FastQC)

Created a program that took an CSV file containing normalized fluorescent intensitites for each nucleotide position in each read as input and used the Viterbi algorithm to compute the most likely sequence of nucleotides for each read in FASTA format. The packages in Python used include numpy, HMM, pandas, matplotlib.  

Also created a program that took a FASTQ file and printed the distribution of quality scores across all reads. To summarize the results, a mean distribution was created and saved as a .png file.

## 2:[Suffix array/tree construction and Binary Search Algorithm with FASTA file(s)](https://github.com/sgallegosperez8/CompGen/tree/main/small-genome-assembly)  

A python script was created to intake a genome as a FASTA file, construct into a a suffix array/tree, and read a pattern file and print all positions where that read occurs in the genome using the bisections algorithm

## 3: [Small Genome Assembly, Alignment and Variant Calling](https://github.com/sgallegosperez8/CompGen/tree/main/suffix-tree_binary-search)
With the reference genome of H5N1, from the SRA. 

QualityQC, Trimming and Small Genome Assembly

Assemble the above genome using the same k-mer size using (a) velvet (https://github.com/dzerbino/velvet Links to an external site.), (b) SPADES (https://github.com/ablab/spades Links to an external site.). Then compare the final assemblies with respect to a reference genome for your species (see full list of species below) using QUAST. Discuss the assemblies with respect to the 3C criteria (completeness, contiguity, and correctness) - once again, note that you will not receive any points if you don't discuss your results (40 points). Submit your QUAST output (as .html file) along with a Word/TeX file with your discussion. Using https://www.ncbi.nlm.nih.gov/sra/SRX32003758[accn] 
