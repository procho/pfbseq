#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson

#quick and dirty poisson using 'trybam2.sam' results from coverage and depth reads
#data = np.loadtxt('trybam2_depth.txt', delimiter='\t')
#data_col3 = data[:,2]
#data = pd.read_csv('trybam2_depth.txt', header=None,sep ='\t')
#print(data.head)
with open('chr1_depth_column.txt', 'r') as col:
    col_list = col.read().split('\n')
    col_list = col_list[:-1]
    col_ints = [int(x) for x in col_list]
#hist = np.histogram(col_ints, bins=75)
x = np.arange(0, 75, 1)
y = poisson.pmf(x, 0.789, 0)
z = poisson.cdf(x, 0.789, 0)
a = poisson.ppf(x, 0.789, 0)
plt.plot(x, y)
plt.plot(x, z)
plt.plot(x, a)
plt.hist(col_ints, bins = 75)
plt.show()

value = poisson.ppf(5, 0.789, 0)
print(value)
