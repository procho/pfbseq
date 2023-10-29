#!/usr/bin/env python3
import sys

# Read in .gaf files (contains full list of human and mouse GO record, downloaded from https://current.geneontology.org/products/pages/downloads.html
# remove less reliable evidence
species = sys.argv[1]

if species == "mouse":
	fh = open("mgi.gaf","r")
elif species == "human":
	fh = open("goa_human.gaf","r")
else:
	print(f"Usage: {sys.argv[0]} <species = \"mouse\" or \"human\">")
	exit(1)

line_count = 1

GOdict = {} # Dictionary, keys are unique genes, values are lists of UNIQUE(!!) GO identifiers associated with genes.

less_reliable_annotations = ["HDA","HGI","HMP","HEP","RCA","IEP"]
Gene_set = set()

for line in fh:
	if line.startswith("!"):
		line_count += 1
		# Skipping header lines
		continue
	else:
		# Extracting gene name and GO from this line
		line.rstrip()
		info_list = line.split("\t")
		Gene = info_list[2]
		Gene_set.add(Gene)

		GO = info_list[4]
		evidence = info_list[6]
		
		if evidence not in less_reliable_annotations:
			
			# Add this GO entry to GOdict. Create Gene key in GOdict if not already exist
			if Gene not in GOdict:
				GOdict[Gene] = set()

			GOdict[Gene].add(GO)
			
		line_count += 1

fh.close()
#print(GOdict)
# Write gene names and GO terms to new file
fh_out = open(f"GO_parser_{species}.out", "w")

#fh_out.write("test 123")
for gene in GOdict:
	fh_out.write(f"{gene}\t")
	for GO in GOdict[gene]:
		fh_out.write(f"{GO}\t")
	fh_out.write("\n")

#fh_out.close()	

