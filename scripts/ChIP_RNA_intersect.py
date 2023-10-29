#!/usr/bin/env python3

import pandas as pd
import sys
from GF_functions import *

# This script will take four inputs: An annotated ChIP-seq file of peaks, a DEseq output file, a log2fc cutoff, and a p value cutoff

peak_file = sys.argv[1] 
DEG_file = sys.argv[2] 
log2FC = sys.argv[3] 
pval = sys.argv[4] 

if len(sys.argv) != 5:
  print('Please include the correct four inputs!')
  sys.exit(1)

deg_set_list = subset_degs(DEG_file, log2FC, pval)

intersect_list=set_intersection(deg_set_list, peak_file)

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
Percent upregulated genes bound by TF: {len(intersect_list[0])/len(deg_set_list[0])}\nPercent downregulated genes bound by TF: {len(intersect_list[1])/len(deg_set_list[1])}''')
  summary_write.write(f'\n_____________________________________________________________________________________________________________________\n')



# create a few different output graphs showing overlap between degs and TF-bound genes!

# First making a venn diagram showing intersections
upregulated_genes_set = deg_set_list[0]
downregulated_genes_set = deg_set_list[1]
chip_df = pd.read_csv(peak_file, sep='\t', header = None)
bound_genes_set = set(chip_df[0])

venn3([upregulated_genes_set, downregulated_genes_set, bound_genes_set], ("Up_Genes","Down_Genes","Bound_Genes"))

plt.savefig('venn_diagram.png')

# Now making a bar graph showing intersections
percent_up_bound = len(intersect_list[0])/len(deg_set_list[0])
percent_down_bound = len(intersect_list[1])/len(deg_set_list[1])

y_list = [percent_up_bound,percent_down_bound]
x_list = ['Up Genes','Down Genes']

fig,ax = plt.subplots(figsize = (8,8))
ax.bar(x_list, y_list, color = 'maroon', edgecolor = 'black')
ax.set_ylabel("% DEGs bound by TF")
ax.grid(b=True, color='grey', linestyle = '-.',linewidth = 0.5, alpha=0.2)
ax.set_title('RNA-seq Differential Genes vs ChIP-seq Peaks')
fig.savefig('bar_graph.png')




