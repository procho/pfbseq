Project Proposal and Presentation Guide
==============

Data processing	Demultiplexing of CITE-seq raw base call files into FASTQ files was performed using 10x Genomics cellranger mkfastq (v.3.0.2).
To perform comparisons between these methodologies, all scRNA-seq and scATAC-seq libraries were downsampled to 200 million FASTQ entries (200M of each sequenced read; ~20,000 reads per expected cell recovered at ~10,000 cell expected recovery), and all ADT libraries were downsampled to 48M reads (Limited by the available read count for 1 experiment; ~4,800 reads per expected cell recovered).
For CITE-seq data, RNA alignment was performed using cellranger (v5.0.0) with the –include-introns parameter against 10x Genomics reference refdata-gex-GRCh38-2020-A. To analyze CITE-seq ADT libraries we used BarCounter to compute the counts of ADT barcodes for each 10x cell barcode
genome_build: Ensembl GRCh38-3.0.0
processed_data_files_format_and_content: Single cell CSVs (Commas separated), and HDF5 file (Hierarchical Data Format).
 	
Project Proposal: 

Overall goal of the project is to analyze CITE-seq data (combination of single cell transcripts and antibody based ADT label of the cells). This is a nice immune cell dataset with a lot of different cell types for analysis. If the project is too small, additional modalities (TEA-seq)(FACS (25 color) can be compiled onto the data as in the paper below.

This project will enable the team to get their hands on this cutting edge data set, explore the python tools that enable analysis, and get data into the correct format to feed into the tools for analysis. Importantly, the publication helps to give direction, but there are likely new ways to do this--and new tools that have been developed in this rapidly changing field.

Inputs:
 The data is in two parts, an scRNAseq dataset (H5 format)--that is used to store multidimensional arrays. The second file is an .csv file that contains the antibody based ADT (barcodes). The goal of the project will be to develop modules that will enable parsing and manipulation of this data, and merging of the datasets to define cell types with ADT components and overlay the scRNAseq data, once the cell types are defined. There is also a Flo-data set. Data is published as part of a larger paper/dataset:Swanson E, Lord C, Reading J, Heubeck AT et al. Simultaneous trimodal single-cell measurement of transcripts, epitopes, and chromatin accessibility using TEA-seq. Elife 2021 Apr 9;10. PMID: 33835024 https://pubmed.ncbi.nlm.nih.gov/33835024/

GEO info: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5123955 CITE-seq, scATAC-seq, single-cell and single-nuclei 10x Multiome ATAC + Gene Expression, and TEA-seq libraries from the same PBMC sample in parallel.

DEMO Examples:
https://cellgeni.github.io/notebooks/html/new-10kPBMC-Scanpy.html
https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html


Data: https://github.com/AllenInstitute/aifi-swanson-teaseq#cite (The .h5 file has the RNA-seq data) The .csv file has the matrix of surface marker ADT counts

Outputs:
Parsed data files that enable analysis
Alignments
Probably most important is the ability to couple the ADT data to improve the Cell Typing in the scRNA seq data.
Data can be analyzed as in the manuscript:
UMAPs; Heat Maps, Other annotation Other resources, cell types fractions of total cells 
***the scientist involved in analysis of this data can provide consultation if needed   


Potential Challenges:

This is likely to be a challenging, but rewarding project. One challenge may be ensuring the data is properly normalized etc. in its current format. Another challenge may be in subsampling the data for the development. Other challenges may be deciding how to implement the outputs. Probably most challenging would be to implement the multimodal components of this dataset- if time allows.  If this is beyond the capabilities of the team, a co-author of the manuscript has offered to advise if needed.


CITE-seq of leukapheresis-purified, 0.01% digitonin permeabilized, FACS neutrophil-depleted PBMCs, method comparison
Identifiers
BioSample: SAMN18105809; GEO: GSM5123955
Organism
Homo sapiens (human)
cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Euarchontoglires; Primates; Haplorrhini; Simiiformes; Catarrhini; Hominoidea; Hominidae; Homininae; Homo
Attributes
source name	peripheral blood mononuclear cells
donor_id	HMN200910
donor_source	BioIVT
purification_type	leukapheresis
preparation_type	whole cells
preparation_conditions	N/A
cleanup_type	FACS
cleanup_conditions	Dead Cell/Debris/Neutrophil Depletion
library_type	CITE-seq
library_type	CITE-seq
Links
GEO Sample GSM5123955
BioProject
PRJNA663623 Simultaneous trimodal single cell measurement of transcripts, epitopes, and chromatin accessibility using TEA-seq
Retrieve all samples from this project

Submission
Allen Institute, Allen Institute For Immunology; 2021-03-02
Accession: SAMN18105809 ID: 18105809
BioProject GEO DataSets

## Strong final project presentations integrate:
- Live demos (on command line or in Jupyter notebook)
- Code walk-throughs/highlights
- Discussion of lessons learned
- Future directions/improvements


## Time line
- Proposals due: Friday the 27th at 2pm. (smr -at- stowers.org prochnik -at- gmail.com  )
- Project sign up: Friday 5PM 
- Halloween Party: Friday 9-12
- Halloween After Party: Friday 12...... 
- Project sign up: Saturday 8:30AM - 9AM   
- First group meeting: Saturday 9AM
- Final Presentaions: Monday 3-5pm
- Banquet: Monday 5:45pm
- Dance party Monday: 9pm

1.  HDF5 files viewed in python.
