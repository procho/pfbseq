#!/usr/bin/env python3

import sys
import pandas as pd
import argparse

##Define function that will take in cluster table with name and pvalues and spit out significant gene names



def split_clusters(filename,outdir = './'): #define name and input for function

    df = pd.read_table(filename) #read file as pandas data frame


    for i in range(0,len(df.columns),2): # loop over list from 0 to number of columns by 2 - [0,2,4,6 ...]
 
        cluster_df = df.iloc[:,[i,i+1]] #define cluster_df as a two column df for the current(i) cluster number 
        cluster_df = cluster_df[cluster_df.iloc[:,1] < 0.0001] # filte for only pvalues smaller than 10e-3
        cluster_df = cluster_df.drop(cluster_df.columns[1],axis = 1) #drop the pvalue column
        out_name = cluster_df.columns[0][0] # grab the cluster number for output file
        cluster_df.to_csv(f'{outdir}/cluster{out_name}_geneNames.tsv',index=False,header=False) # save 1 column df without index and without header

def main():

    progname = sys.argv[0]
    
    parser = argparse.ArgumentParser(description="")
    
    parser.add_argument('file', help = '.tsv file that needs to be split and filtered based on p-values')    
    parser.add_argument("-o", "--outdir",metavar = 'outdir', required = False, default = './', help = 'optional: path to desired output directory')
    args = parser.parse_args()    
    
    split_clusters(args.file, args.outdir)

    
    print('Complete!\nFiles Generated: clusterX_geneNames.tsv')
    sys.exit(0)

if __name__ == '__main__':
    main()



