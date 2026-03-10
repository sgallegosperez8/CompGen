import Bio
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from hmmlearn import hmm
from Bio import SeqIO


#Question 1: Create HMM Base Caller

#open file to store sequences
seq_reads = open("assignment1.fasta","w")

#loading csv data
data = pd.read_csv("assignment1.csv")

#filtered df to only focus on the '1' positions of the reads to use for read lengths
index_pos =data[data["Position"]== 1]

#stored read lengths in an array
index_positions = data.index.get_indexer_for(index_pos.index)
idx_p = index_positions[1:]

nuc = ['A','T','G','C']
observations = data[["Intensity_A","Intensity_T","Intensity_G","Intensity_C"]].values
trans_matrix = np.array([
    [0.9,0.033,0.033,0.033],
    [0.033,0.9,0.033,0.033],
    [0.033,0.033,0.9,0.033],
    [0.033,0.033,0.033,0.9]
])

#normalized the transition matrix
trans_matrix = trans_matrix/trans_matrix.sum(axis=1,keepdims=True)

#model
model = hmm.GaussianHMM(n_components=4,covariance_type="diag")

#parameters
model.startprob_ = np.array([0.25,0.25,0.25,0.25])
model.transmat_ = trans_matrix

model.init_params = 'mc'  # Only initialize means and covariances
model.params = 'mc'

#fit the model
model.fit(observations)
log_prob, state_seq = model.decode(observations, algorithm="viterbi")
predicted_path = [nuc[state] for state in state_seq]

#use the idx_p(read lengths) to split the reads by indices and into a file called 'assignment1.fasta'
for i in range(0,len(idx_p)):
    if i == 0:
        first_seq =''.join(predicted_path[:(idx_p[i]-1)])
        seq_reads.write(f'>seq {i+1}\n{first_seq}\n')
    elif i < 23:
        middle_seq =''.join(predicted_path[idx_p[i-1]:idx_p[i]-1])
        seq_reads.write(f'\n>seq {i+1}\n{middle_seq}\n')
    elif i ==23:
        last_2_seq =''.join(predicted_path[idx_p[i-1]:idx_p[i]-1])
        seq_reads.write(f'\n>seq {i+1}\n{last_2_seq}\n')
        last_seq =''.join(predicted_path[idx_p[i]:])
        seq_reads.write(f'\n>seq {i+2}\n{last_seq}\n')
        
#file close
seq_reads.close()




#Question 2: Creating a mean distribution based off the fastq file

#parse through fastq file and change ASCII character to quality score 
quality_char = []
for record in SeqIO.parse("assignment1.fastq","fastq"):
    quality_char.append(record.letter_annotations["phred_quality"])

#use a for loop to place all quality score for each char associated in a dictionary then dataframe
data_dict = {}
for i in range(len(quality_char)):
    data_dict[f"seq{i+1}"] = quality_char[i]

data_df = pd.DataFrame(data_dict)

#use for loop to retrieve all values for "each" row in the data frame and retrieve the mean 
sum_list = []
for i in range(len(data_df["seq1"])):
    sum_list.append(data_df[i:i+1].values.mean())

# create an x axis and the means will be the y axis 
col_names = np.arange(1,101)
sum_list = np.array(sum_list)

#plot means on a bar graph
plt.bar(col_names,sum_list)

plt.title('Mean Quality Score Across Sequence Read Position')
plt.xlabel('Position in read (bp)')
plt.ylabel('Quality Score')

plt.show()
