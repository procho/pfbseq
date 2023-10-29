#!/usr/bin/env python3

import pandas as pd
import sys

# This script will take a multi-sheet input excel file and will export a tab separated text file

# Inputs:
####### 1- excel file
####### 2- sheet name in excel file
####### 3- Desired output file name (without extension)

if len(sys.argv) != 4:
  print('Usage: tsv_from_xlsx.py <excel file> <sheet name> <output name>')
  sys.exit(1)


excel_bed = sys.argv[1]
sheet_name = sys.argv[2]
output_file = sys.argv[3]

df_bed = pd.read_excel(excel_bed, sheet_name = sheet_name)


df_bed.to_csv(f'{output_file}.txt', sep='\t', index=False, header=False)

