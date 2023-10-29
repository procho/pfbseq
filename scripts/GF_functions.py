#!/usr/bin/env python3

import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn3


# This function will take a tab separated file of DEGs (gene	log2fc	pval) and will return two files
# of up or downregulated genes based on specified log2fc and pvalue cutoffs as well as a list of sets
# where list[0] is significantly upregulated genes and list[1] is downregulated genes

def subset_degs(deg_file, log2fc, pval):
  # Read in the DEseq2 table to pandas
  deg_df = pd.read_csv(deg_file, sep='\t', header=None, names=['gene','log2fc','pval'])
  
  all_genes = deg_df['gene']
  all_genes.to_csv('all_genes.txt', sep='\t', index=False, header=False)

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



# This function will take a list of sets of DEGs where list[0] is upregulated genes and list[1] is
# downregulated genes. It will also take a file containing annotated chip seq peaks where the gene
# names are in column[0]. 
# The script will return a list of sets where list[0] is upregulated bound genes and list[1] is
# downregulated bound genes 

def set_intersection(set_list, anno_chip):
  chip_df = pd.read_csv(anno_chip, sep='\t', header = None)
  peak_set = set(chip_df[0])
  
  intersect_list = []
  
  for DE_set in set_list:
    intersect = DE_set.intersection(peak_set)
    intersect_list.append(intersect)

  return intersect_list
    


# This function will take three sets and three names and will create a three-way venn diagram
# showing the intersection between all three sets with the appropriate name labels

def venn_diagram(set1, set2, set3, name1, name2, name3):
  venn3([set1,set2,set3], (name1, name2, name3))
  plt.savefig('venn_diagram.png')


# This function will take two lists, one with X values and one with Y values, and will make a bar graph

def bar_graph(x_list, y_list):
  fig,ax = plt.subplots(figsize = (8,8))
  ax.bar(x_list, y_list, color = 'maroon', edgecolor = 'black')
  ax.set_ylabel("% DEGs bound by TF")
  ax.grid(b='visible', color = 'grey', linestyle = '-.', linewidth = 0.3, alpha = 0.2)
  ax.set_title("RNA-seq Differential Genes vs ChIP-seq Peaks")
  fig.savefig('bar_graph.png')  
