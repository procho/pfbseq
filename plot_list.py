#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as py
from scipy.stats import poisson
values = []

with open('trybam2_depth.txt', 'r') as depth:
    for line in depth:
        values.append(line.split('\t')[2])
    print('done')
print(values)
plt.plot(values)
plt.ylabel('read depth')
plt.show()
plt.savefig('rd_dpth_chr1.png')
