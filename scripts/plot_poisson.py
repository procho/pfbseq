#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as py
from scipy.stats import poisson
values = []
file = sys.argv[1]
base_depth = sys.argv[2]
total_coverage = 
d_file = f'{file}.depth.txt'
c_file = f'{file}.cov.txt'

#prep list of depth values for regular graph
with open(d_file, 'r') as depth:
    for line in depth:
        values.append(line.split('\t')[2])

#get values to run poisson
with open(c_file, 'r') as cover:
    for line in cover:
        
print(values)
plt.plot(values)
plt.ylabel('read depth')
plt.show()
plt.savefig('rd_dpth_chr1.png')


~                                                                                            
~                 
