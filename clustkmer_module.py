#!/usr/bin/env python
# coding: utf-8



# Import libraries 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import numba
#import py3Dmol
import pandas as pd 
import random
import math 

from tqdm.notebook import tqdm
import sklearn
import copy 
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score




# Define a function to read a .tsv file and return its contents as a numpy array
def read_pssm(fp):
    with open(fp) as file:
        lines = file.readlines()
    return np.array([[int(score) for score in line.rstrip('\n').split('\t')[2:]] for line in lines[2:]])



## Function to read pssm file and msa file 

def process_files(pssm_path, msa_path):
    
    # Call the function with the first file path and store the result in the variable pssm 
    pssm_ss = read_pssm(pssm_path)
            
    # Open the MSA file
    aln_seqs = []
    with open(msa_path) as file:
       # Iterate over each line in the file
       for line in file:
           # Skip any comment lines
           if line.startswith(('#', '//')):
               continue
           # Split the line into two parts at the first space character
           _, seq = line.split(' ')
           # Remove any trailing newline characters from the sequence and append it to the list
           aln_seqs.append(seq.rstrip('\n'))
            
    return pssm_ss, aln_seqs 





## Function to define the score matrix based on pssm and msa file 
alphabet = 'ACDEFGHIKLMNPQRSTVWY'

def kmer_scores(pssm_path, msa_path):
    # Call process_files to read and process the files
    pssm, aln_seqs = process_files(pssm_path, msa_path)

    # Initialize a 2D numpy array filled with zeros. The dimensions of the array are determined 
    # by the length of the PSSM and the number of sequences.
    scores = np.zeros((len(pssm), len(aln_seqs)))

    # Iterate over each position in the PSSM
    for i in range(2, len(pssm)-2):
        # Iterate over each sequence in aln_seqs
        for j, seq in enumerate(aln_seqs):
            # Initialize a score variable to zero
            score = 0
            # Iterate over five positions relative to the current position in the PSSM
            for k in [-2, -1, 0, 1, 2]:
                # Get the amino acid at the current position in the sequence
                letter = seq[i + k]
                # If the amino acid is a gap, add np.nan to the score. Otherwise, add the PSSM value for the 
                # amino acid at the corresponding position in the PSSM.
                if letter == '-':
                    score += np.nan
                else:
                    score += pssm[i + k, alphabet.index(letter)]
            # Assign the score to the corresponding position in the scores array
            scores[i, j] = score
    # Return the scores array
    return scores






## Function to show the profile plot based on the scores extracted from pssm and msa file 
def profile_show(pssm_path, msa_path):
    scores_ss = kmer_scores(pssm_path, msa_path)
    im = plt.imshow(scores_ss.T, vmin=-35, vmax=45)
    plt.title('Profiles scores of 3Di 5-mers')
    plt.ylabel('MSA sequence')
    plt.xlabel('MSA position')
    plt.colorbar(im)
    plt.tight_layout()
    plt.show()
    return im 


## Function to receive neg/pos samples and make them numerical 
def process_fasta(file_path):
    
    fasta = []
    with open(file_path) as fasta_file:
        cnt = 0
        for line in tqdm(fasta_file, disable=True):
            if not line.startswith('>'): 
                fasta.append(np.array([alphabet.index(l) for l in line.rstrip('\n')]))
                cnt +=1
        print("Number of cluster members:", cnt)
    return fasta

## Function to extract kmer properties including kmer position, kmer counts and score threshold for all seqs 
def kmer_properties(scores_ss, fasta,pssm_ss):
    # For each column, take the 10 percentile score as threshold
    score_threshold = np.percentile(np.nan_to_num(scores_ss, -100), 10, axis=1)
    # Remove columns with a threshold of 0
    kmer_pos = np.where(score_threshold)[0]
    score_threshold = score_threshold[kmer_pos]
    # Count for each columns the hits above the threshold
    kmer_thr_counts = np.zeros_like(kmer_pos)
    
    for seq in tqdm(fasta, disable=True):
        count_kmers5(pssm_ss, kmer_pos, score_threshold, seq, kmer_thr_counts)
    
    return score_threshold, kmer_pos, kmer_thr_counts 


## Function to count occurance of kmers in sequences (dataset)
@numba.njit(parallel=False)
def count_kmers5(pssm, kmer_pos_lst, score_threshold_lst, seq, counts):
    already_counted = np.zeros(len(kmer_pos_lst))
    for seq_pos in numba.prange(2, len(seq)-2):        
        for kmer_i in range(len(kmer_pos_lst)):
            score = 0
            for k in [-2, -1, 0, 1, 2]:
                score += pssm[kmer_pos_lst[kmer_i] + k, seq[seq_pos + k]]
            if score >= score_threshold_lst[kmer_i] and already_counted[kmer_i] != 1:
                counts[kmer_i] += 1
                already_counted[kmer_i] = 1



## Function to convert ouputs of kmer_properties to dataframe 
def kmer_hit_df(kmer_thr_counts, kmer_pos, score_threshold, aln_seqs):
    # Initialize an empty list to store the information
    info_list = []

    idx = np.argsort(kmer_thr_counts)

    for k, (pssm_pos, cnt, thr) in enumerate(zip(kmer_pos[idx], kmer_thr_counts[idx],
                                                score_threshold[idx])):
        # Save the information in a dictionary
        info_dict = {
            '3Di': aln_seqs[0][pssm_pos-2:pssm_pos+3],
            #'AA': aln_aa_seqs[0][pssm_pos-2:pssm_pos+3],
            'Hits': cnt,
            'Threshold': thr,
            'Position': pssm_pos,
            'idx' : idx[k]
        }
        # Append the dictionary to the list
        info_list.append(info_dict)

    # Convert the list of dictionaries into a DataFrame
    hits_df = pd.DataFrame(info_list)
    
    return hits_df




## Write a function that receives the +/- samples and select some portion of them 

def select_and_label(neg_sample, pos_sample, neg_prc, pos_prc):
    
    # Randomly select 'num_values1' elements from 'list1'
    neg_num = round(neg_prc*len(neg_sample))
    neg_selected = random.sample(neg_sample, neg_num)
    
    # Assign label 0 to these elements 
    neg_labels = [0]*len(neg_selected)
    
    # Randomly select 'num_values2' elements from 'list2'
    pos_num = round(pos_prc*len(pos_sample))
    pos_selected = random.sample(pos_sample, pos_num)
    
    # Assign label 1 to all elements in 'selected_list2'
    pos_labels = [1]*len(pos_selected)
    
    # Combine two lists 
    combined_list = neg_selected + pos_selected
    combined_labels = neg_labels + pos_labels
    
    # Pair each element with its respective label
    # combined_list = list(zip(combined_list, combined_labels))
    combined_list = pd.DataFrame({
        'Element': combined_list,
        'Label': combined_labels 
    })
    return combined_list 



## Function for classification of each sequence based on number of kmers (n), number of matched kmers (threshold)
## and distance between the kmers 
def best_kmer_set2(n,thresh,distance, hits_df, pssm_ss, combined_list,scores_ss):
    
#     # make new list of kmers with distance conditioning 
#     dist_idx = [0]
#     pos_val = [hits_df['Position'][0]]
    
    
#     for i in range(1,len(hits_df)):
#         pos = hits_df['Position'][i]                       # Take the position of kmer
#         pos_subtract = [abs(j-pos) for j in pos_val]       # Subtract the pos of current kmer from the position of 
#                                                            # all of previous kmers         
#         if all([k>=distance for k in pos_subtract]):       # If the subtraction result for all of the values is more than 
#             pos_val.append(pos)                            # than *distance* add the current kmer to the list 
#             dist_idx.append(i)
    
#     seqs_subset_df_dist = hits_df.iloc[dist_idx].reset_index()
#     print(hits_df)                        ######### ADDED FOR DEBUGGING ############# 
#     print(seqs_subset_df_dist)            ######### ADDED FOR DEBUGGING ############# 
#     if n>len(seqs_subset_df_dist):
#         print("The legth of kmer set size,",n,"exceeds the size of the available kmers,",len(seqs_subset_df_dist))
#         return None
    prc = 100000
    l = 0
    while l == 0 and prc<=1000000:
        seqs_ten_prc = hits_df[hits_df['Hits']<= prc]    # Select those kmers with frequecy less than prc
        l = len(seqs_ten_prc) 
        if l==0:
            prc +=100000
        if l==1:
            thresh = 1

    if l < n:
        n=l
    seqs_subset_df_dist = seqs_ten_prc
    
    # Initialize list to store labels 
    labels=np.zeros(len(combined_list), dtype=int)
    
    # Get kmer positions from pssms_ss
    score_threshold = np.percentile(np.nan_to_num(scores_ss, -100), 10, axis=1)
    kmer_pos = np.where(score_threshold)[0]
    score_threshold = score_threshold[kmer_pos]
    
    # Count for each columns the hits above treshold 
    kmer_thr_counts = np.zeros_like(kmer_pos)

    # Extract positions 
    pos = seqs_subset_df_dist['Position'][:n]

    score_threshold = seqs_subset_df_dist['Threshold'][:n]
    
    # Loop through each element in the combined_list 
    for idx, seq in enumerate(tqdm(combined_list['Element'][:], disable=True)):
        # Initialize counts 
        counts = np.zeros_like(kmer_pos)

        # Count hits for kmer
        count_kmers5(pssm_ss, np.array(pos), np.array(score_threshold), seq, counts)

        # Assign 1 if more than/equal threshold kmers are present, otherwise 0
        hits_one = counts.sum()
        
        if hits_one>= thresh:
            labels[idx]=1
                       
    
    # Return updated combined_list 
    return labels



## Function to calculate the metrics based on the predicted_labels 
def metric_calc(pred_list, heatmap=False):
    # Assuming y_true is your true labels and y_pred is your predicted labels
    y_true = pred_list['Label']
    y_pred = pred_list['Pred_Label']


    # Calculate confusion matrix
    #print('y_true:',y_true)                         ######### ADDED FOR DEBUGGING ############# 
    #print('y_pred:',y_pred)                         ######### ADDED FOR DEBUGGING ############# 
    cm = confusion_matrix(y_true, y_pred)
    #print(cm)
    tn, fp, fn, tp = cm.ravel()
    
    if heatmap:
        # Create a heatmap with labels and values
        sns.heatmap(cm, annot=True, fmt='', cmap='Blues')
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.show()

    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)
    specificity = tn/(tn+fp)
    g_mean = math.sqrt(recall*specificity)
    f1 = f1_score(y_true, y_pred)


#     print(f'Accuracy: {accuracy}')
#     print(f'precision: {precision}')
#     print(f'Sensitivity (TP rate): {recall}')
#     print(f'Specificity (TN rate): {specificity}')
#     print(f'Geometric mean of Sens-Spec: {g_mean}')
#     print(f'f1: {f1:.3f}')
    
    return accuracy, precision, recall, specificity, g_mean, f1
    




