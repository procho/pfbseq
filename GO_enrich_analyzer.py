#!/usr/bin/env python3

import sys
from scipy.stats import fisher_exact
from GO_counter import *

def fold_enrich(sample_file,ctrl_file, speciesDict):
	
	# get sample GO counts, store as a dictionary 
	GOdict_sample = GO_count(sample_file, speciesDict)
	
	
	## Get number of total genes in sample
	fh_sample = open(sample_file, "r")
	gene_num_sample = 0
	for line in fh_sample:
		line = line.rstrip()
		if line:
			gene_num_sample += 1
	
	fh_sample.close()

	
	# get control GO counts, store as a dictionary
	GOdict_ctrl = GO_count(ctrl_file,speciesDict)
	
	## Get number of total genes in control	
	fh_ctrl = open(ctrl_file, "r")
	gene_num_ctrl = 0
	for line in fh_ctrl:
		line = line.rstrip()
		if line:
			gene_num_ctrl += 1
	
	fh_ctrl.close()
		
	# For each GO, calculate fold enrichment, and statistical test
	
	GO_results_Dict = {}
	GO_list = []

	for GO in GOdict_ctrl:  ## Assuming sample gene list is a subset of control gene list
		
		GO_list.append(GO)
		GO_results_Dict[GO] = {}

		if GO in GOdict_sample:
			GOcount_sample = GOdict_sample[GO]
		else:
			GOcount_sample = 0
		
		GOcount_ctrl = GOdict_ctrl[GO]

		fold_enrichment = (GOcount_sample/gene_num_sample) / (GOcount_ctrl/gene_num_ctrl)
		
		# Fisher's exact test, obtain p value
		# simon: Jan 2024
		# fisher exact test is expecting counts of the same thing, not a mixture of GO annots
		# and gene counts. Need to take a look at the lines below
		print('Something is not quite right here. Check the code.')
		exit(1)
		a = GOcount_sample
		b = GOcount_ctrl
		c = gene_num_sample - GOcount_sample
		d = gene_num_ctrl - GOcount_ctrl
		table = [[a,b],[c,d]]
		Fisher_result = fisher_exact(table)
		p = Fisher_result.pvalue
		
		GO_results_Dict[GO]["fold_enrich"] = fold_enrichment
		GO_results_Dict[GO]["pvalue"] = p

	GO_list_sorted = sorted(GO_list, key=lambda x:GO_results_Dict[x]["pvalue"])
	GO_results_Dict_sorted = {}
	
	for go in GO_list_sorted:
		GO_results_Dict_sorted[go] = GO_results_Dict[go]	

	return(GO_results_Dict_sorted)



def main():
	pass 
	# Get species Gene:GO dictionary
	#speciesDict = GeneGoDict(f"GO_parser_human.out")
	#foldchange_human = fold_enrich("fake_sample_human.txt", "fake_ctrl_human.txt", "GO:0005515",speciesDict)
	#foldchange_mouse = fold_enrich(fake_sample_mouse, fake_ctrl_mouse, "??")
	#print(foldchange_human)

#main()
