#!/usr/bin/env python3
import sys
from subprocess import call
import os

### Write functions for each step of the workflow:
# bam file ---> samtools depth output
# depth output -(discard lines below read threshhold)-> depths above threshhold file
# thresh file ---> list of peaks ---> bed file of peaks
print(sys.argv[1], sys.argv[2])

def bam_to_depth(bam_file): #input bam file, runs samtools depth
    call(f"samtools depth -o {bam_file[:-4] + '.depth'} {bam_file}", shell=True)


def depth_to_thresh(depth_file): #input depth output, runs awk. awk makes new file with
#all lines that contain read depth values greater than 5.
    print('nothing!')
    cmd = "awk '{if($3>5)print$0}'" + f" {depth_file} > {depth_file[:-6] + '.thresh.depth'}"
    os.system(cmd)

def thresh_to_peaks(thresh_file):
    num_lines = os.system(f"wc -l {thresh_file}")
    peaks = []
    chr_x = ''
    start = 0
    end = 0
    x = 0 #assign max depth var
    line_num = 0
    with open(thresh_file, 'r') as t_file:
        for line in t_file: #current script LOSES THE LAST PEAK because addition of peaks
#is only triggered by moving to a new peak.
            line_num += 1
            line = line.split('\t')
            pos = int(line[1])
            chr_x = line[0]
            base_depth = int(line[2])
            if base_depth > x: #keep largest base read depth in each peak
                x = base_depth
            if line_num == num_lines:
                peaks.append([chr_x, str(start), str(end), str(x)])
            elif pos > end+1:
                peaks.append([chr_x, str(start), str(end), str(x)])
                start = pos
                end = pos
                x = 0
            elif pos == end+1:
                end += 1
            elif pos <= end and chr_x == line[0]:
                print('This file is not properly ordered. It needs to be ordered by chromosome and chromosome location.')
                exit(1)
    out_file = f'{thresh_file[:-12]}bed'
    with open(out_file, 'w') as out:
        for i_list in peaks:
            line = "\t".join(i_list) + '\n'
            out.write(line)
    
funcs = {'bam_to_depth': bam_to_depth, 'depth_to_thresh': depth_to_thresh, 'thresh_to_peaks' : thresh_to_peaks}

#def main():

if __name__ == '__main__':
    funcs[sys.argv[1]](sys.argv[2])
