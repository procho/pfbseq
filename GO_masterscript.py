#!/usr/bin/env python3
import sys
import os
from GO_enrich_analyzer import *

species = sys.argv[1]

# Get species GO dictionary
speciesDict = GeneGoDict(f"GO_parser_{species}.out")

if species == "mouse":

	# For mouse RNAseq/ChIPseq data
	ctrl_file = "./Danni/all_genes.txt"
	upsample_file = "./Danni/upregulated_bound_genes.txt"
	downsample_file = "./Danni/downregulated_bound_genes.txt"

	#upregulated genes
	GO_results_up = fold_enrich(upsample_file, ctrl_file, speciesDict)

	up_fh = open("GO_mouse_RNAseq_upregulated.out","w")
	up_fh.write("GO_term\tfold_change\tpvalue\n")

	for GO in GO_results_up:
		fold_change = GO_results_up[GO]["fold_enrich"]
		pvalue = GO_results_up[GO]["pvalue"]
		up_fh.write(f"{GO}\t{fold_change}\t{pvalue}\n")

	up_fh.close()

	#downregulated genes

	GO_results_down = fold_enrich(downsample_file, ctrl_file, speciesDict)

	down_fh = open("GO_mouse_RNAseq_downregulated.out","w")
	down_fh.write("GO_term\tfold_change\tpvalue\n")

	for GO in GO_results_down:
		fold_change = GO_results_down[GO]["fold_enrich"]
		pvalue = GO_results_down[GO]["pvalue"]
		down_fh.write(f"{GO}\t{fold_change}\t{pvalue}\n")

	down_fh.close()		

elif species == "human":
	
	# For human PBMC ADT data
	ctrl_file = "split_all_PMBC_genes.csv"

	for i in range(0,19):
		sample_file = f"./clusterGenes/cluster{i}_geneNames.tsv"
		GO_results = fold_enrich(sample_file,ctrl_file,speciesDict)
		
		fh_w = open(f"GO_PBMC_cluster{i}_geneNames.out","w")
		fh_w.write("GO_term\tfold_change\tpvalue\n")
	
		for GO in GO_results:
			fold_change = GO_results[GO]["fold_enrich"]
			pvalue = GO_results[GO]["pvalue"]
			fh_w.write(f"{GO}\t{fold_change}\t{pvalue}\n")
		
		fh_w.close()

else:
	# Help message
	print(f"Usage: {sys.argv[0]} <species = \"human\" or \"mouse\">")		
	exit(1)

