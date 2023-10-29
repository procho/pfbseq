#!/usr/bin/env python3

import sys
import subprocess
from goatools import obo_parser

def GO_id2term(go_id):

	go_obo = "/Users/pfb09/final_project/pfbseq/data/go-basic.obo"
	go = obo_parser.GODag(go_obo)
	go_term = go[go_id]
	return go_term

def GeneGoDict(filename):
	# A function that generate dictionary with genes and associated GO terms. IN: GO_parser.out file(human or mouse)
	GO_Dict = {}
	fh = open(filename,"r")
	for line in fh:
		line = line.rstrip()
		splitline = line.split("\t")
		gene = splitline[0]
		GO_list = splitline[1:]
		GO_Dict[gene] = GO_list
	
	return(GO_Dict)
	

def GO_count(filename, speciesDict):

	# Read gene names from file
	Genelist_fh = open(filename, "r")
	GO_list_full = []
	
	# Generate a list of all GOs associated with gene indicated in input file. Stored in GO_list_full.
	for gene in Genelist_fh:
		gene = gene.rstrip()
		
		if gene not in speciesDict:
			continue

		GO_list = speciesDict[gene]
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
	GO_dict = GeneGoDict("GO_parser_human.out")
	GO_count_dict = GO_count("fake_ctrl_human.txt",GO_dict)
	print(GO_count_dict)
#main()
	
