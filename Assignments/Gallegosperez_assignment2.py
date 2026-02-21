from Bio import SeqIO




for record in SeqIO.parse("Assignment2_refgenome.fasta", "fasta"):
    print (record.id)

    










