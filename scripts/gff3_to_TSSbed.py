#!/usr/bin/env python3

import re
import sys


# This function will take an input GFF3 file, and will locate the start
# of each transcript (taking strand into consideration)

# It will export a .bed formatted file with the start position of each 
# identified transcription start site!

# gff3 file downloaded from gencode
def gff3_to_tssbed(GFF3_file):
  with open(GFF3_file, 'r') as GFF_read, open('TSS.bed', 'w') as GFF_write:
    for line in GFF_read:
      line = line.rstrip()
      #Get rid of header lines
      if not line.startswith('#'):
        #Make a list from each line, split on tab
        line_list = line.split('\t')
        #Extract only transcript entries from line_list[2]
        if line_list[2] == 'transcript':
          label_string = line_list[8] 
          #Extract gene name from long string in line_list[8]
          match = re.search(r'gene_name=([^;]+);', label_string)
        
          # Setting variables to write export file
          chrom = line_list[0]
          gene_sym = match.group(1) 
          strand = line_list[6] 
          if line_list[6] == '+':
            # For + strand, take the start position and account for bed-gff3 file differences
            bed_start = int(line_list[3]) - 1
            bed_end = int(line_list[3])
          else:
            # For - strand, take the end position and account for bed-gff3 file differences
            bed_start = int(line_list[4]) - 1
            bed_end = int(line_list[4])
   
          GFF_write.write(f'{chrom}\t{bed_start}\t{bed_end}\t{gene_sym}\t{strand}\n')     


