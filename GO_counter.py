#!/usr/bin/env python3

import sys
from goatools import obo_parser

def GO_id2term():

	go_obo = "go-basic.obo"
	go = obo_parser.GODag(go_obo)
	return go

def GeneGoDict(filename):
	# A function that generate dictionary with genes and associated GO terms. IN: GO_parser.out file(human or mouse)
	id2term_dict = GO_id2term()
	
	GO_Dict = {}
	fh = open(filename,"r")
	for line in fh:
		line = line.rstrip()
		splitline = line.split("\t")
		gene = splitline[0]
		GO_list = splitline[1:]
		GO_termlist = [id2term_dict[GOid] for GOid in GO_list]
		GO_Dict[gene] = GO_termlist
	
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
	
