#!/bin/zsh


# This script takes two args - the chip bed file and the TSS bed file
chip_bed=$(basename $1 .bed)
TSS_bed=$(basename $2 .bed)

# Sort these bed files based on chromosome and start position 
sort -k1,1 -k2,2n $1 > ${chip_bed}.sorted.bed
sort -k1,1 -k2,2n $2 > ${TSS_bed}.sorted.bed

# Run bedtools closest and pull out the 7th column using awk which should have closest annotated genes
bedtools closest -a ${chip_bed}.sorted.bed -b ${TSS_bed}.sorted.bed | awk '{ print $7 }' > annotated_bound_genes.txt

rm ${chip_bed}.sorted.bed ${TSS_bed}.sorted.bed
