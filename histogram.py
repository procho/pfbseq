#!/usr/bin/env python3

from scipy.stats import poisson
import matplotlib.pyplot as plt
import numpy
import sys

# For plotting a histogram of a samtools depth file.

file_name = sys.argv[1]
y_axis = sys.argv[2] #Name of the y-axis

if len(sys.argv) > 3:
    x_axis = sys.argv[3] #If 3rd arg added, x_axis will be named, and actual chromosome vals will be tracked
else:
    x_axis = None

d_list = []
loc_list = []

with open(file_name, 'r') as file:
    for line in file:  
        line = line.split('\t')
        d_list.append(line[2])
        if bool(x_axis) == True:
            loc_list.append(line[1])

#if bool(loc_list) == True:
#    plt.plot(loc_list, d_list)
#    plt.ylabel(y_axis)
#    plt.xlabel(x_axis)
#    plt.show()
#else:
#    plt.plot(d_list)
#    plt.ylabel(y_axis)
#    plt.show()
#plt.savefig(f'x_y_plot_{file_name[:-3]}')
#plt.hist(loc_list)
#n, bins, patches = patches = plt.hist(x, 75, density = True, facecolor = 'g', alpha=0.75)

