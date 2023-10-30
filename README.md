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
(description of peak calling algorithm + script)

* **Peak Annotation:**
ChIP-seq peaks are annotated based on a user input GFF3 genome annotation file.
First, the GFF3 file is parsed using ```gff3_to_TSSbed.py``` which extracts
transcription start site information and gene names from the GFF3 file. Then,
ChIP-seq peaks are annotated to the closest transcription start site using
bedtools in ```annotate_bed_to_TSS.sh```. 
 
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





# Part 2: 
