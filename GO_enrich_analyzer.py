#!/usr/bin/env python3

import sys
from scipy.stats import fisher_exact
from GO_counter import *

def fold_enrich(sample_file,ctrl_file, GO, species):
	
	# get sample GO count, and number of total genes
	GOdict_sample = GO_count(sample_file, species)
	GOcount_sample = GOdict_sample[GO]
	
	fh_sample = open(sample_file, "r")
	gene_num_sample = 0
	for line in fh_sample:
		line = line.rstrip()
		if line:
			gene_num_sample += 1
	
	fh_sample.close()
	print(GOcount_sample)	
	print(gene_num_sample)
	
	# get control GO count, and number of total genes
	GOdict_ctrl = GO_count(ctrl_file,species)
	GOcount_ctrl = GOdict_ctrl[GO]
	
	fh_ctrl = open(ctrl_file, "r")
	gene_num_ctrl = 0
	for line in fh_ctrl:
		line = line.rstrip()
		if line:
			gene_num_ctrl += 1
	
	fh_ctrl.close()
	print(GOcount_ctrl)
	print(gene_num_ctrl)
		
	#print(gene_num_sample, gene_num_ctrl)
	# calculate fold enrichment, return this ratio
	fold_enrichment = (GOcount_sample/gene_num_sample) / (GOcount_ctrl/gene_num_ctrl)
	
	# Fisher's exact test, obtain p value
	a = GOcount_sample
	b = GOcount_ctrl
	c = gene_num_sample - GOcount_sample
	d = gene_num_ctrl - GOcount_ctrl
	table = [[a,b],[c,d]]
	Fisher_result = fisher_exact(table)
	p = Fisher_result.pvalue
	
	return(fold_enrichment,p)



def main():
	foldchange_human = fold_enrich("fake_sample_human.txt", "fake_ctrl_human.txt", "GO:0005515","human")
	#foldchange_mouse = fold_enrich(fake_sample_mouse, fake_ctrl_mouse, "??")
	print(foldchange_human)

main()
