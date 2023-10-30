#!/usr/bin/env python3

import pandas as pd
import sys
from GF_functions import *
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import subprocess
from gff3_to_TSSbed import *
import argparse
from argparse import RawTextHelpFormatter
from bam_to_peaks import *
import os

###### Formatting and documenting input arguments #######

parser = argparse.ArgumentParser(description=f"This is a wrapper script for comparison of RNAseq and ChIPseq data.\n\nIf a bam or a sam file is provided from a ChIP-seq experiment, a custom peak calling algorithm will be run. Otherwise a pre-made bed file may be entered with called ChIP-seq peaks.\n\nAn RNAseq file should be provided containing differential expression data, it should be a tab-delimited file containing columns for gene IDs, fold change (log2 transformed), and p values. If desired, the user may designate custom pval and log2fc cutoffs for downstream analysis. Otherwise, reasonable default values will be utilized.\n\nThis script will annotate called peaks (either from an input bed file or generated from a bam file) based on the distance to the nearest gene TSS. It will then output many files and graphs representing overlap and intersections between significant differentially expressed genes (either up or downregulated) to called ChIP-seq peaks. Finally, the script will output enriched GO terms for either upregulated bound genes or downregulated bound genes.", formatter_class = RawTextHelpFormatter)

parser.add_argument("-G", required = True, help='Input GFF3 annotation file for genome of ChIP/RNAseq data', dest = 'GFF3_file')

parser.add_argument("-r", required = True, help='Input differential gene expression file. Required column format = gene_ID  Log2FC  p-value', dest = 'RNAseq_file')

parser.add_argument("-b", required = False, help='Optional: supply input ChIP seq bed file, if this is not provided, user must provide a .bam file for peak calling', dest = 'bed_file')

parser.add_argument("-B", required = False, help='Optional: supply input bam file for peak calling before intersection. If this is not provided, user must provide a ChIPseq .bed file', dest = 'bam_file')

parser.add_argument("-f", required = False, help='Optional: supply log2fc cutoff for differentially expressed gene list. Default: 0.585', dest = 'log2fc')

parser.add_argument("-p", required = False, help='Optional: supply p value cutoff for differentially expressed gene list. Default: 0.05', dest = 'pval')

args = parser.parse_args()


# Quit and report error if user inputs both a bed file and a bam file
if args.bed_file and args.bam_file:
  print('Please only input a bed file OR a bam file!')
  sys.exit(1)

# Quit and report error if user inputs neither a bed file or a bam file
if not args.bed_file and not args.bam_file:
  print('Must input either a bed file or a bam file!')
  sys.exit(1)


######### Assigning input args to variables ########

GFF3_file = args.GFF3_file
DEG_file = args.RNAseq_file

if args.pval:
  pval = args.pval
else:
  pval = 0.05

if args.log2fc:
  log2FC = args.log2fc
else:
  log2FC = 0.585

if args.bam_file:
  bam_file = args.bam_file

if args.bed_file:
  peak_file = args.bed_file


# If bam file is input, run peak caller!
if args.bam_file:
  bam_to_depth(bam_file)
  depth_to_thresh(f"{bam_file[:-4] + '.depth'}")
  thresh_to_peaks(f"{bam_file[:-4] + '.thresh.depth'}")
  peak_file = f"{bam_file[:-4] + '.bed'}"  


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


# Run GO analysis on upregulated bound genes and downregulated bound genes!
cmd2 = f'./GO_masterscript.py mouse'
cmd2_run = subprocess.run(cmd2, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

if cmd_run.returncode != 0:
  print(f'GO analysis step failed')
  exit(2)
