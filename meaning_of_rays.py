# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:15:33 2016

@author: cave
"""

import numpy
import scipy
import scipy.stats
#http://stats.stackexchange.com/questions/99777/calculating-the-variance-of-the-histogram-of-a-grayscale-image
data = numpy.load("psfs.npy")
current_slice = 0
while current_slice < 20:
    under_review = data[current_slice,:]
    x_values = numpy.linspace(-30,30,4096)
    pre_mean = under_review*x_values
    mean=numpy.sum(x_values*under_review)/(1000*1000)
    varience=numpy.sum(((x_values-mean)*(x_values-mean))*under_review)/(1000*1000)
    range_data = scipy.stats.norm(mean, varience).interval(0.80)
    print range_data
    print current_slice*100+100
    current_slice = current_slice+1