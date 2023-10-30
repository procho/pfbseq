## Single Cell RNA-seq to UMAP walkthrough

### Overview

In this overview we will highlight the various functions of scanpy in the analysis of single-Cell RNA-seq data and the creation of UMAP graphs.

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

* Regressing out data

```python
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
```

* Clustering

```python
sc.pp.neighbors(adata, n_neighbors=5, n_pcs=7)
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

* Leiden Clustering

Resolution will change number of clusters found:

```python
sc.tl.leiden(adata, resolution = 1.2, key_added = "leiden_1.2")
sc.tl.leiden(adata, resolution = 1.5, key_added = "leiden_1.5")
sc.tl.leiden(adata, resolution = 2.0, key_added = "leiden_2.0")
```
```
running Leiden clustering
    finished: found 14 clusters and added
    'leiden_1.2', the cluster labels (adata.obs, categorical) (0:00:00)
running Leiden clustering
    finished: found 15 clusters and added
    'leiden_1.5', the cluster labels (adata.obs, categorical) (0:00:00)
running Leiden clustering
    finished: found 19 clusters and added
    'leiden_2.0', the cluster labels (adata.obs, categorical) (0:00:00)
```

```python
sc.pl.umap(adata, color=['leiden_1.2', 'leiden_1.5', 'leiden_2.0'], legend_loc='on data', wspace = 0.25, legend_fontsize=10)
```

```python
sc.pl.umap(adata, color=['leiden', 'CD4', 'CD8A', 'TCF7L2'])
```

* Can do clustering based on kmeans or Agglomerative Clustering

```python
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
X_pca = adata.obsm['X_pca'] 

kmeans = KMeans(n_clusters=15, random_state=0).fit(X_pca) 
adata.obs['kmeans15'] = kmeans.labels_.astype(str)
kmeans = KMeans(n_clusters=20, random_state=0).fit(X_pca) 
adata.obs['kmeans20'] = kmeans.labels_.astype(str)
kmeans = KMeans(n_clusters=25, random_state=0).fit(X_pca) 
adata.obs['kmeans25'] = kmeans.labels_.astype(str)

sc.pl.umap(adata, color=['kmeans15', 'kmeans20', 'kmeans25'], legend_loc='on data', wspace = 0.25, legend_fontsize=10)
```
```python
from sklearn.cluster import AgglomerativeClustering

cluster = AgglomerativeClustering(n_clusters=15, affinity='euclidean', linkage='ward')
adata.obs['hclust15'] = cluster.fit_predict(X_pca).astype(str)

cluster = AgglomerativeClustering(n_clusters=20, affinity='euclidean', linkage='ward')
adata.obs['hclust20'] = cluster.fit_predict(X_pca).astype(str)

cluster = AgglomerativeClustering(n_clusters=25, affinity='euclidean', linkage='ward')
adata.obs['hclust25'] = cluster.fit_predict(X_pca).astype(str)

sc.pl.umap(adata, color=['hclust15', 'hclust20', 'hclust25'], legend_loc='on data', wspace = 0.25, legend_fontsize=10)
```




### Antibody walkthrough

* Reading in .csv file with cell and antibody counts

```python

```

