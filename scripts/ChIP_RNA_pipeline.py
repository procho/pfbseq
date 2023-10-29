#!/usr/bin/env python3

import pandas as pd
import sys
from GF_functions import *
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import subprocess
from gff3_to_TSSbed import *

# This script will take five inputs: A GFF3 file to use for annotation, A ChIP-seq bed file of peaks, a DEseq output file, a log2fc cutoff, and a p value cutoff
GFF3_file = sys.argv[1]
peak_file = sys.argv[2] 
DEG_file = sys.argv[3] 
log2FC = sys.argv[4] 
pval = sys.argv[5] 

if len(sys.argv) != 6:
  print('Please include the correct four inputs!')
  sys.exit(1)

# Run function to create TSS bed file from GFF3 file
gff3_to_tssbed(GFF3_file)

# Use subprocess to run bash script to export annotated peaks
cmd = f'./annotate_bed_to_TSS.sh {peak_file} TSS.bed'
cmd_run = subprocess.run(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

if cmd_run.returncode != 0:
  print(f'Bedtools sort and annotation step failed')
  exit(2)


# Run function to subset input DEG file
deg_set_list = subset_degs(DEG_file, log2FC, pval)

# Intersect annotated peaks and DEGs
intersect_list=set_intersection(deg_set_list, 'annotated_bound_genes.txt')

# create output text files containing up or downregulated genes bound by TF
with open('upregulated_bound_genes.txt', 'w') as up_write:
  for bound_gene in intersect_list[0]:
    up_write.write(f'{bound_gene}\n')
  
with open('downregulated_bound_genes.txt','w') as down_write:
  for bound_gene in intersect_list[1]:
    down_write.write(f'{bound_gene}\n')


with open('intersect_summary.txt', 'w') as summary_write:
  summary_write.write(f'_____________________________________________________________________________________________________________________\n')
  summary_write.write(f'''Number of upregulated genes: {len(deg_set_list[0])}\nNumber of downregulated genes: {len(deg_set_list[1])}\n\n
Number of upregulated genes bound by TF: {len(intersect_list[0])}\nNumber of downregulated genes bound by TF: {len(intersect_list[1])}\n\n
Fraction upregulated genes bound by TF: {len(intersect_list[0])/len(deg_set_list[0])}\nFraction downregulated genes bound by TF: {len(intersect_list[1])/len(deg_set_list[1])}''')
  summary_write.write(f'\n_____________________________________________________________________________________________________________________\n')



# create a few different output graphs showing overlap between degs and TF-bound genes!

# First making a venn diagram showing intersections
upregulated_genes_set = deg_set_list[0]
downregulated_genes_set = deg_set_list[1]
chip_df = pd.read_csv('annotated_bound_genes.txt', sep='\t', header = None)
bound_genes_set = set(chip_df[0])

venn_diagram(upregulated_genes_set, downregulated_genes_set, bound_genes_set, 'Up_Genes','Down_Genes','Bound_Genes')



# Now making a bar graph showing intersections
percent_up_bound = len(intersect_list[0])/len(deg_set_list[0])
percent_down_bound = len(intersect_list[1])/len(deg_set_list[1])

y_list = [percent_up_bound,percent_down_bound]
x_list = ['Up Genes','Down Genes']

bar_graph(x_list, y_list)

