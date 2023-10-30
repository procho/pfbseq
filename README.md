# Multi-Omics Toolkit
A toolkit for the intersection and analysis of omics datasets



## Part 1: Integration of ChIP-seq and RNA-seq datasets

### Description and Rationale:
As omics experiments grow more and more commonplace, the need for tools to integrate these various data modalities has grown. This wrapper script allows for quick and easy integration of two common datasets, bulk RNA-seq and ChIP-seq.


### Overview of Pipeline: ```ChIP_RNA_pipeline.py```
<img src='./images/biorender_workflow.png' width='700' height='300'/> 

```
usage: ChIP_RNA_pipeline.py [-h] -G GFF3_FILE -r RNASEQ_FILE [-b BED_FILE] [-B BAM_FILE] [-f LOG2FC]
                             [-p PVAL]

This is a wrapper script for comparison of RNAseq and ChIPseq data.

If a bam or a sam file is provided from a ChIP-seq experiment, a custom peak calling algorithm will be run.
Otherwise a pre-made bed file may be entered with called ChIP-seq peaks.

An RNAseq file should be provided containing differential expression data
and should be a tab-delimited file containing columns for gene IDs, fold change (log2 transformed), and p values.
If desired, the user may designate custom pval and log2fc cutoffs for downstream analysis.
Otherwise, reasonable default values will be utilized.

This script will annotate called peaks (either from an input bed file or generated from a bam file)
based on the distance to the nearest gene TSS. It will then output many files and graphs
representing overlap and intersections between significant differentially
expressed genes (either up or downregulated) to called ChIP-seq peaks.
Finally, the script will output enriched GO terms for either upregulated bound genes or downregulated bound genes.

options:
  -h, --help      show this help message and exit
  -G GFF3_FILE    Input GFF3 annotation file for genome of ChIP/RNAseq data
  -r RNASEQ_FILE  Input differential gene expression file. Required column format = gene_ID  Log2FC  p-value
  -b BED_FILE     Optional: supply input ChIP seq bed file, if this is not provided, user must provide a .bam file for peak calling
  -B BAM_FILE     Optional: supply input bam file for peak calling before intersection. If this is not provided, user must provide a ChIPseq .bed file
  -f LOG2FC       Optional: supply log2fc cutoff for differentially expressed gene list. Default: 0.585
  -p PVAL         Optional: supply p value cutoff for differentially expressed gene list. Default: 0.05
```


* **Peak Calling:**
Input sam or bam files are run through ```samtools depth```, a function which
returns the locations and read depth of each covered base in the input file.
This output is filtered by removing all bases with a read depth below 6, and
consecutive bases are grouped into peaks and output to a bed file.

* **Peak Annotation:**
ChIP-seq peaks are annotated based on a user input GFF3 genome annotation file.
First, the GFF3 file is parsed using ```gff3_to_TSSbed.py``` which extracts
transcription start site information and gene names from the GFF3 file. Then,
ChIP-seq peaks are annotated to the closest transcription start site using
[bedtools closest](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html) in ```annotate_bed_to_TSS.sh```. 
 
* **Integration of RNA-seq and ChIP-seq:**
The input RNAseq gene list is parsed to identify up and downregulated genes
based on user-defined cutoffs. These lists are then integrated with annotated
ChIP-seq peaks. A number of summary statistics are reported in a summary.txt
file, as well as lists of significant DE genes, and lists of overlap between
upregulated or downregulated genes and annotated peaks. A bar graph and a pie
chart are also generated to visualize overlap between these datasets:

<img src='./images/venn_diagram.png' width='400' height='300'/>        <img src='./images/bar_graph.png' width='300' height='300'/>


* **GO Pathway Analysis of DEGs and ChIP-seq Peaks**
Gene Ontology pathway analysis is performed on lists of overlapping annotated
peaks with either upregulated or downregulated genes. Results are exported
in a table format with relevant statistics:

```
GO_term fold_change     pvalue
GO:0070161      level-03        depth-03        anchoring junction [cellular_component] 4.821307290568866       3.3283419434173545e-09
GO:0005886      level-03        depth-03        plasma membrane [cellular_component]    1.8713643260863588      1.7589741581267755e-08
GO:0016020      level-02        depth-02        membrane [cellular_component]   1.565372257152236       3.452880176773256e-07
GO:0015629      level-06        depth-06        actin cytoskeleton [cellular_component] 4.968408197939896       7.903768061919691e-07
GO:0042995      level-02        depth-02        cell projection [cellular_component]    2.5492854843103747      1.2241457808974728e-06
GO:0009897      level-03        depth-03        external side of plasma membrane [cellular_component]   3.988385412504329       5.073967906570723e-06
```

The environment used to run these scripts can be downloaded using the included [env.yaml](https://github.com/procho/pfbseq/blob/main/share/RNA_ChIP.yaml) file: 
```
conda env create --file RNA_ChIP.yaml --name RNA_ChIP
```
Note: Samtools must be installed separately!





## Single Cell RNA-seq to UMAP walkthrough

### Overview

In this overview we will highlight the various functions of scanpy in the analysis of single-Cell RNA-seq data and the creation of UMAP graphs for the scRNA-seq data and Antibody counts.

### RNA walkthrough

We will be using the numpy, pandas, scanpy, matplotlib packages all downloadable through conda. Environment set up instructions can be found in jupyterSetUp.md

* Importing Packages into notebook

```python
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
```
* Reading in the HD5 file with the data

```python
adata = sc.read_10x_h5('/Users/pfb16/final_project/pfbseq/scSEQ_data/GSM5123955_X066-RP0C1W1_leukopak_perm-cells_cite_200M_rna_counts.h5')#replace with path to your dataset
adata.var_names_make_unique()
```
* Annotating mitochondrial genes

```python
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```
* Filtering data sets based on desired parameters

```python
sc.pp.filter_cells(adata, min_genes=200)# filters out cells with few transcripts
sc.pp.filter_genes(adata, min_cells=3) # filters out genes expressed in few cells
adata = adata[adata.obs.n_genes_by_counts < 2500, :] #Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

* Normalizing data and accounting for highly variable data

```python
sc.pp.normalize_total(adata, target_sum=1e4) #Logarithmize the data

sc.pp.log1p(adata) #Identify highly-variable genes.

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) # extracting highly variable genes
adata.var.highly_variable.value_counts()

adata = adata[:, adata.var.highly_variable] #Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.

```

* Make a file with all gene names that will be clustered

```python
out_file = 'all_PBMC_genelist.txt'
adata.var[['gene_ids']].to_csv('./all_PBMC_genelist.csv')
```

* Regressing out mitochondrial data

```python
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
```

* Clustering

```python
sc.pp.neighbors(adata, n_neighbors=5, n_pcs=7)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.8)
```
* Generating UMAP

```python
with plt.rc_context({'figure.figsize': (8, 8)}):
    sc.pl.umap(adata, color=['leiden'],legend_loc='on data')
```

![image](https://github.com/procho/pfbseq/assets/110238030/563b6db8-2127-4961-9eca-8176d449d350)





### Antibody walkthrough

* Reading in .csv file with cell and antibody counts

```python
protein = sc.read_csv('./share/GSM5123955_X066-RP0C1W1_leukopak_perm-cells_cite_48M_adt_counts_fixed.csv')
protein.var_names_make_unique()
```

* Graph of antibody counts

```python
sc.pl.highest_expr_genes(protein, n_top=40, )
```
![image](https://github.com/procho/pfbseq/assets/110238030/0e08b7b5-ed3a-449f-9c49-a19b9846b223)

* Normalization and clustering

```python3
sc.pp.log1p(protein)
sc.pp.pca(protein, n_comps=20)
sc.pp.neighbors(protein, n_neighbors=30)  
sc.tl.leiden(protein,resolution = 0.8)
```

* Creating UMAP

```python
sc.tl.umap(protein)
with plt.rc_context({'figure.figsize': (8, 8)}):
    sc.pl.umap(protein, color=['leiden'],legend_loc='on data')
```

![image](https://github.com/procho/pfbseq/assets/110238030/7c1e2e6b-0ded-420d-8700-147e37e96d55)

