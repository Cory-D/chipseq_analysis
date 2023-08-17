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

#loading files, receive parameters, and providing help
ap = argparse.ArgumentParser()
ap.add_argument('-i','--input',required=True,type=str,help="Input: file in CSV format\n")
ap.add_argument('-r','--rows',required=False,type=int,default=0,help="The number of rows you would like to capture (first row ignored).")
ap.add_argument('-p','--percent',required=False,type=int,default=0,help="The percent of rows you would like to capture (first row ignored).")
ap.add_argument('-t','--times',required=False,type=int,default=1,help="The number of times to sample the starting file.")

args = vars(ap.parse_args())

file_input = args['input']
rows_to_choose = args['rows']
percent_choice = args['percent']
sample_times = args['times']

#load input file

input_table_DF = pd.read_csv(file_input,sep = ',',header=0) # load comma delimited text as dataframe
list_of_column_headers = list(input_table_DF.columns.values.tolist())
total_rows = len(input_table_DF.index)

#remove CSV from file input

file_input = file_input[:-4]

#decide on approach

if percent_choice != 0 and rows_to_choose == 0:
    count = 0
    fraction_choice = percent_choice/100
    for i in range(sample_times):
        table_sample = input_table_DF.sample(frac = fraction_choice)
        print(table_sample)
        save_filename = file_input + '_random_' + str(i) + '.csv'
        table_sample.to_csv(save_filename,sep=',',index=False,header=True)

elif percent_choice == 0 and rows_to_choose != 0:
    count = 0
    column_count = len(list_of_column_headers)
    for i in range(sample_times):
        input_table_DF = input_table_DF.sample(frac=1) #shuffle the spreadsheet
        table_sample = input_table_DF.iloc[0:rows_to_choose,0:column_count]
        print(table_sample)
        table_sample.columns = list(list_of_column_headers)
        save_filename = file_input + '_random_' + str(i) + '.csv'
        table_sample.to_csv(save_filename,sep=',',index=False,header=True)
else:
    print("Neither number of rows, nor percentage of rows, has been selected as a parameter.")

