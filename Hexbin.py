#######################################################################################################################
#       A script to generate a density plot of intron coverage values in one techniques versus another 
#
#                   To run, place in folder with comparative data values as two column .csv :
#
#                                             nohup Hexbin.py & 
#
#     Coded by Robert Beagrie, small modifications by Martin Larke - University of Oxford 30-Aug-2019
#######################################################################################################################

#!/usr/bin/env python

import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import numpy as np
 
# read mean RPKM values into a dataframe
reviews = pd.read_csv('1_hexbin.txt', sep='\t')

# Set figure size
plt.figure(figsize=(8,6))
plt.hexbin(x=reviews.total, y=reviews.su,

           gridsize=60,        # this parameter sets the number of hexagons (higher integeter means more hexagons across each axis)
           xscale='log', yscale='log', cmap=plt.cm.BuPu)


# Log scaled hexbins don't set the axis limits correctly for some reason

# auto scaled
plt.xlim(reviews.total.min(), reviews.total.max())
plt.ylim(reviews.su.min(), reviews.su.max())

# custom scaled e.g xlim(left,right)
#plt.xlim(0.001, 1000)
#plt.ylim(0.001, 1000)

plt.xlabel('Total RNA')
plt.ylabel('PAM')

# Make a colorbar with sensible ticks
cb = plt.colorbar()
cb.ax.set_ylabel('Number of genes', rotation=90)

# print some statistics
pearson_r, p_val = stats.pearsonr(np.log10(reviews.total),

                                  np.log10(reviews.su))
print 'pearson r'
print pearson_r
print 'p_val'
print p_val
 
plt.savefig('my_hexbin.pdf')
