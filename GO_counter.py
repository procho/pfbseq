#!/usr/bin/env python3

import sys
import subprocess
from goatools import obo_parser

def GO_id2term(go_id):

	go_obo = "/Users/pfb09/final_project/pfbseq/data/go-basic.obo"
	go = obo_parser.GODag(go_obo)
	go_term = go[go_id]
	return go_term

def GO_count(filename, species):

	if species == "mouse" or species == "human":
		GO_out_filename = f"GO_parser_{species}.out"
	else:
		print(f"Usage: {sys.argv[0]} <GeneList_filename> <species>")
		exit(2)
	
	# Read gene names from file
	Genelist_fh = open(filename, "r")
	GO_list_full = []
	
	# Generate a list of all GOs associated with gene indicated in input file. Stored in GO_list_full.
	for gene in Genelist_fh:
		gene = gene.rstrip()
		record = subprocess.check_output(f"grep -w {gene} {GO_out_filename}",shell=True)
		record = record.decode("utf-8")
		record = record.rstrip()
		GO_list = record.split("\t")[1:]
		GO_list_full.extend(GO_list)

	# Count each GO in the full GO list
	
	GO_count_dict = dict()
	for GO in GO_list_full:
		if GO not in GO_count_dict:
			GO_count_dict[GO] = 1
		else:
			GO_count_dict[GO] += 1
	return(GO_count_dict)

def main():
	GO_count("fake_ctrl_human.txt", "human")
	#result = GO_count(sys.argv[1],sys.argv[2])
	#print(result)
	#print(GO_id2term("GO:0048527"))

main()
	
