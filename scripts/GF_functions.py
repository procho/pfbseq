#!/usr/bin/env python3

import pandas as pd
import sys

### This function will take three inputs:
### 1) A DEseq file with DEGs, log2FC, pval, 2) desired log2fc cutoff, 3) desired p value cutoff

### The function will return:
### A list of sets where list[0] = upregulated genes and list[1] = downregulated genes according to the user input cutoffs

# The script will also return a summary .txt. file with the number and identities of differentially up- and downregulated genes according to input cutoffs 





def subset_degs(deg_file, log2fc, pval):
  # Read in the DEseq2 table to pandas
  # Must be in this format! with the column names in this order....
  deg_df = pd.read_csv(deg_file, sep='\t', header=None, names=['gene','log2fc','pval'])

  upregulated = deg_df[(deg_df['log2fc'] > float(log2fc)) & (deg_df['pval'] < float(pval))]
  upregulated.to_csv('upregulated_genes.txt',sep='\t')



  downregulated = deg_df[(deg_df['log2fc'] < -(float(log2fc))) & (deg_df['pval'] < float(pval))]
  downregulated.to_csv('downregulated_genes.txt',sep='\t')  

  up_set = set(upregulated['gene'])
  down_set = set(downregulated['gene'])
  
  set_list = []
  set_list.append(up_set)
  set_list.append(down_set)
  
  return set_list
 

def set_intersection(set_list, anno_chip):
  # Have a list of upregulated and downregulated genes
  # Also include annotated ChIP seq peaks (with genes in column 1)
  chip_df = pd.read_csv(anno_chip, sep='\t', header = None)
  peak_set = set(chip_df[0])
  intersect_list = []
  for DE_set in set_list:
    intersect = DE_set.intersection(peak_set)
    intersect_list.append(intersect)

  return intersect_list
    
 

def main():
  input_peak_file = sys.argv[1]
  input_deg_file = sys.argv[2]
  input_log2fc = sys.argv[3]
  input_pval = sys.argv[4]
  
  deg_set_list = subset_degs(input_deg_file, input_log2fc, input_pval)
  
  print(f'Number of upregulated genes: {len(deg_set_list[0])}\nNumber of downregulated genes: {len(deg_set_list[1])}')

  intersect_list=set_intersection(deg_set_list, input_peak_file)
  print(intersect_list)

if __name__ == '__main__':
  main()
  
