#!/usr/bin/env python
# coding: utf-8
# Author: Cory Dunn
# Institution: University of Helsinki
# License: None
# Author Email: cory.dunn@helsinki.fi

#import numpy and pandas and argparse

import pandas as pd #import pandas
import numpy as np #import numpy
import argparse

##prepare data for analysis

#loading files, receive parameters, and providing help
ap = argparse.ArgumentParser()
ap.add_argument('-g','--genomic_sites',required=True,type=str,help='The list of genomic sites should be within a CSV file with precisely these headers in the first row:\n'\
    'Chromosome | Start_Nucleotide | End_Nucleotide | GeneID | TSS\n')
ap.add_argument('-c','--chipseq_hit_sites',required=True,type=str,help='The list of ChipSeq hits should be within a CSV file with precisely these headers in the first row:\n'\
    'Chromosome | Start_Nucleotide | End_Nucleotide | PeakID')
ap.add_argument('-w','--window_distance',required=False,type=int,default=1000,help='The window size (bp) used to scan near the transcriptional start site (default: 1000 bp).')
ap.add_argument('-j','--jump_distance',required=False,type=int,default=500,help='The step movement of the window when scanning near the transcriptional start site (default: 500 bp).')
ap.add_argument('-r','--range_distance',required=False,type=int,default=10500,help='The furthest distance from the TSS to the right or left border of the region of interest (TSS -range-> | (default: 10500).')

args = vars(ap.parse_args())
print('Running on genomic_sites file: {}'.format(args['genomic_sites']))
print('Running on chipseq_hit_sites file: {}'.format(args['chipseq_hit_sites']))
genomic_file = args['genomic_sites']
chipseq_file = args['chipseq_hit_sites']
window_size = args['window_distance']
jump_size = args['jump_distance']
rang_around_TSS = args['range_distance']

#load list of selected genomic sites for comparison
random_gene_full_table = pd.read_csv(genomic_file,sep = ',') # load comma delimited text as dataframe
random_gene_full_table.dropna # remove any rows with any NA
print(random_gene_full_table)

#convert columns of genomic dataframe to individual lists
random_gene_chr_list = random_gene_full_table['Chromosome'].tolist()
random_gene_start_list = random_gene_full_table['Start_Nucleotide'].tolist()
random_gene_end_list = random_gene_full_table['End_Nucleotide'].tolist()
random_gene_id_list = random_gene_full_table['GeneID'].tolist()
random_gene_TSS_list = random_gene_full_table['TSS'].tolist()

#load list of ChipSeq hits

chipseq_hit_full_table = pd.read_csv(chipseq_file,sep = ',') # load comma delimited text as dataframe

print(chipseq_hit_full_table)

#convert columns of chipseq dataframe to individual lists

chipseq_hit_chr_list = chipseq_hit_full_table['Chromosome'].tolist()
chipseq_hit_start_list = chipseq_hit_full_table['Start_Nucleotide'].tolist()
chipseq_hit_end_list = chipseq_hit_full_table['End_Nucleotide'].tolist()
chipseq_hit_peakID_list = chipseq_hit_full_table['Peak_ID'].tolist()

##compare genomic and chipseq lists

#initialize lists
list_of_ChipSeq_hits_NO_OVERLAP_in_random = []
list_of_ChipSeq_hits_YES_OVERLAP_in_random = []

#for each chipseq hit, test all entries in genomic list
for i, chiphit in enumerate(chipseq_hit_peakID_list):
    for j, randhit in enumerate(random_gene_id_list):
        if chipseq_hit_chr_list[i] == random_gene_chr_list[j] and chipseq_hit_start_list[i] <= random_gene_end_list[j] and chipseq_hit_start_list[i] >= random_gene_start_list[j]:
            list_of_ChipSeq_hits_YES_OVERLAP_in_random.append([chipseq_hit_peakID_list[i],chipseq_hit_chr_list[i],chipseq_hit_start_list[i],chipseq_hit_end_list[i],random_gene_id_list[j],random_gene_chr_list[j],random_gene_start_list[j],random_gene_end_list[j],random_gene_TSS_list[j]])
        elif chipseq_hit_chr_list[i] == random_gene_chr_list[j] and random_gene_start_list[j] >= chipseq_hit_start_list[i] and random_gene_start_list[j] <= chipseq_hit_end_list[i]:
            list_of_ChipSeq_hits_YES_OVERLAP_in_random.append([chipseq_hit_peakID_list[i],chipseq_hit_chr_list[i],chipseq_hit_start_list[i],chipseq_hit_end_list[i],random_gene_id_list[j],random_gene_chr_list[j],random_gene_start_list[j],random_gene_end_list[j],random_gene_TSS_list[j]])
        elif chipseq_hit_chr_list[i] == random_gene_chr_list[j] and random_gene_start_list[j] == chipseq_hit_start_list[i] and random_gene_end_list[j] == chipseq_hit_end_list[i]:
            list_of_ChipSeq_hits_YES_OVERLAP_in_random.append([chipseq_hit_peakID_list[i],chipseq_hit_chr_list[i],chipseq_hit_start_list[i],chipseq_hit_end_list[i],random_gene_id_list[j],random_gene_chr_list[j],random_gene_start_list[j],random_gene_end_list[j],random_gene_TSS_list[j]])
        else:
            list_of_ChipSeq_hits_NO_OVERLAP_in_random.append([chipseq_hit_peakID_list[i],chipseq_hit_chr_list[i],chipseq_hit_start_list[i],chipseq_hit_end_list[i],random_gene_id_list[j],random_gene_chr_list[j],random_gene_start_list[j],random_gene_end_list[j],random_gene_TSS_list[j]])
            

#generate dataframes
NO_OVERLAP = pd.DataFrame(list_of_ChipSeq_hits_NO_OVERLAP_in_random)
YES_OVERLAP = pd.DataFrame(list_of_ChipSeq_hits_YES_OVERLAP_in_random)
NO_OVERLAP.columns = ['chipseq_hit_peakID','chipseq_hit_chr','chipseq_hit_start','chipseq_hit_end','random_gene_id','random_gene_chr','random_gene_start','random_gene_end','random_gene_TSS']
YES_OVERLAP.columns = ['chipseq_hit_peakID','chipseq_hit_chr','chipseq_hit_start','chipseq_hit_end','random_gene_id','random_gene_chr','random_gene_start','random_gene_end','random_gene_TSS']

#find relative position of each hit to TSS

#initialize lists
list_of_hit_midpoints = []
list_of_relative_TSS_distances = []

#midpoints of chipseq hit

count_of_yes = len(list_of_ChipSeq_hits_YES_OVERLAP_in_random)
for i in range(count_of_yes):
    midpoint_chipseq_hit = (YES_OVERLAP.iloc[i]['chipseq_hit_start'] + YES_OVERLAP.iloc[i]['chipseq_hit_end'])/2
    list_of_hit_midpoints.append(int(midpoint_chipseq_hit))

#relative distances of hits from TSS

for i in range(count_of_yes):
    hit_distance = list_of_hit_midpoints[i] - YES_OVERLAP.iloc[i]['random_gene_TSS']
    list_of_relative_TSS_distances.append(hit_distance)

#add to YES_OVERLAP DataFrame

YES_OVERLAP['chipseq_hit_midpoint_distance_from_TSS'] = np.array(list_of_relative_TSS_distances)

#print both and save list of overlaps

print(NO_OVERLAP)
print(YES_OVERLAP)
chipseq_file = chipseq_file[:-4] #remove .csv from filename
chipseq_hit_file_out = 'ChipSeq_hits_of_' + chipseq_file + "_in_" + genomic_file
YES_OVERLAP.to_csv(chipseq_hit_file_out)

#sliding window collection of hits on each size of TSS

window_left = (-1 * rang_around_TSS)
window_right = window_left + window_size
window_list = []
total_sum_across_windows = 0

while window_right <= rang_around_TSS:
    sum_in_window = len([i for i in list_of_relative_TSS_distances if i <= window_right and i >= window_left])
    total_sum_across_windows += sum_in_window
    window_list.append([window_left,window_right,sum_in_window])
    window_left = window_left + jump_size
    window_right = window_right + jump_size

#print and save sum of hits within each window

WINDOW_LIST_DF = pd.DataFrame(window_list)
WINDOW_LIST_DF.columns = ['window_left_of_TSS','window_right_of_TSS','chipseq_hit']
print(WINDOW_LIST_DF)
chipseq_window_file_out = 'Distance_of_ChipSeq_hits_of_' + chipseq_file + "_near_TSS_of_" + genomic_file

WINDOW_LIST_DF.to_csv(chipseq_window_file_out)

#wrap up

print('\nTask completed.\n')
print('Genomic sites used as input: ' + str(len(random_gene_id_list)) + '\n')
print('ChipSeq hit sites used as input : ' + str(len(chipseq_hit_peakID_list)) + '\n')

print('Total number of ChipSeq set matches within genomic site set: ' + str(count_of_yes) + '\n')
print("File output for list of ChipSeq hits in genomic set: '" + chipseq_hit_file_out + "'\n")
print('***Note: Window size, window jump, and range does not apply to the above analysis.\n')

print('Parameters for sliding window across TSS (window width): ' + str(window_size) + '| (window jump): ' + str(jump_size) + '| (range on each side of TSS): ' + str(rang_around_TSS) + "'\n\n")
print("File output for window bins near TSS: '" + chipseq_window_file_out + "'\n")
