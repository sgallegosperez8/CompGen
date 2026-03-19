# Computational Genomic Projects 

## 1: [HMM BaseCaller and FastQC](https://github.com/sgallegosperez8/CompGen/tree/main/HMM-BaseCaller-and-FastQC)

Created a program that took an CSV file containing normalized fluorescent intensitites for each nucleotide position in each read as input and used the Viterbi algorithm to compute the most likely sequence of nucleotides for each read in FASTA format. The packages in Python used include numpy, HMM, pandas, matplotlib.  

Also created a program that took a FASTQ file and printed the distribution of quality scores across all reads. To summarize the results, a mean distribution was created and saved as a .png file.

## 2: [Suffix array/tree construction and Binary Search Algorithm with FASTA file(s)](https://github.com/sgallegosperez8/CompGen/tree/main/small-genome-assembly)  

A python script was created to intake a genome as a FASTA file, construct into a a suffix array/tree, and read a pattern file and print all positions where that read occurs in the genome using the bisections algorithm

## 3: [Small Genome Assembly](https://github.com/sgallegosperez8/CompGen/tree/main/small-genome-assembly)

Based on the accession SRX32003758 (H5N1), raw sequencing reads of a virus were assembled and evaluated.

Below is a summary of the required steps.

1. ***Genome Assembly Process***

Performed de novo assembly using two different algorithms but keeping the k-mer size (=25)constant to ensure a fair comparison.

Velvet: A classic de Bruijn graph assembler. It consists of two main programs: velveth (to hash the reads and create the dataset) and velvetg (to build the graph and produce the assembly).

SPAdes: A more modern, versatile assembler designed for single-cell and standard bacterial datasets. It often performs better with varying coverage and provides built-in error correction.


2. ***Evaluation using QUAST***

With my FASTA files (e.g., contigs.fa from Velvet and contigs_SPAdes.fasta from SPAdes), QUAST (Quality Assessment Tool for Genome Assemblies) was used to compare the two tools. In order to do, a reference genome for the species associated with SRX32003758 was used to generate a comprehensive report.


## 4: [Alignment and Variant Calling](https://github.com/sgallegosperez8/CompGen/tree/main/variant-calling)

A comparative genomic analysis of SARS-CoV-2 was conducted by selecting five distinct samples from the NCBI Virus database, spanning various geographical locations and time points from 2020 through 2026. Using the original Wuhan strain (NC_045512.2) as the primary reference genome, I performed variant calling across the five samples to identify single nucleotide polymorphisms (SNPs) and structural variations. 







