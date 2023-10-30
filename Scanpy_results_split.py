#!/usr/bin/env python3

# To split the output files from Scanpy (Gene ID, ENS ID)

import pandas as pd

data = pd.read_csv('all_PBMC_genelist.csv')

data.drop('gene_ids', inplace=True, axis=1)

print(data.iloc[:,0].to_string(index=False))

data.to_csv('split_all_PMBC_genes.csv',index=False)
