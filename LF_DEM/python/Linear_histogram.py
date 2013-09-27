#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import math
import numpy as np
import string


class LinearHistogram:

    def __init__(self, r_bn, min_r, max_r):
        self.r_bin_nb=r_bn

        self.r_max=max_r
        self.r_min=min_r

        self.r_bsize=(self.r_max-self.r_min)/self.r_bin_nb

        self.histogram = np.zeros(self.r_bin_nb)


    def update(self, distances, value, act='add'):
        
        r_bin = np.array((distances-self.r_min)/self.r_bsize, dtype=int)
        inrange_bin = np.logical_and(r_bin<self.r_bin_nb, r_bin>=0)

        r_bin = r_bin[ inrange_bin ]
        values = np.array(value)[ inrange_bin ]

        if act == 'add':
            unique_indices = np.unique(r_bin)
            aggregated_values = (np.bincount(r_bin, weights = values))[unique_indices]
            self.histogram[unique_indices] += aggregated_values
        elif act == 'replace':
            self.histogram[unique_indices] = values
        else:
            sys.stderr.write("  LinearHistogram :: update(pos, value, coord='sph',act='add') :\n       Unknown act %s".str(act))
            sys.exit(1)


    def normalize(self, norm_factor):

        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            rlo=self.r_bsize*i+self.r_min
            dvol=rlo*self.r_bsize
            self.histogram[i]/=dvol*norm_factor

    
    def print_to(self, stream):
        headline=str(self.r_bin_nb)+' '+str(self.r_max)+'\n'
        stream.write(headline)

        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            out_string=str(r)+" "+str(self.histogram[i])+'\n'
            stream.write(out_string)


    
