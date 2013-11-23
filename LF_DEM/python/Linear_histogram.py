#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import math
import numpy as np
import string


class LinearHistogram:

    def __init__(self, params):

        if params[0] == 'linear':
            [ r_bn, min_r, max_r ] = params[1:]
            self.r_bins = np.linspace(self.min_r, self.max_r, self.r_bn)
        elif params[0] == 'custom':
            self.r_bins = np.array(params[1])

        self.r_bin_nb = self.r_bins.shape[0]-1
        self.r_max = self.r_bins.max()
        self.r_min = self.r_bins.min()
        self.r_bsize=(self.r_max-self.r_min)/self.r_bin_nb

        self.histogram = np.zeros(self.r_bin_nb)


    def update(self, distances, value, act='add'):
        inrange = np.logical_and(distances>self.r_min, distances<self.r_max)
        inrange = np.nonzero(inrange)

        distances = distances [ inrange ]
        values = np.array(value)[ inrange ]

        r_bin = np.digitize(distances,self.r_bins)-1

        if act == 'add':
            unique_indices = np.unique(r_bin)
            aggregated_values = (np.bincount(r_bin, weights = values))[unique_indices]
            self.histogram[unique_indices] += aggregated_values
        elif act == 'replace':
            self.histogram[unique_indices] = values
        else:
            sys.stderr.write("  LinearHistogram :: update(pos, value, coord='sph',act='add') :\n       Unknown act %s".str(act))
            sys.exit(1)


    def normalize(self, norm_factor, norm_type='1d'):

        for i in range(self.r_bin_nb):
            if norm_type == '1d':
                dvol = self.r_bins[i+1] - self.r_bins[i]
            if norm_type == '2d':
                dr = self.r_bins[i+1] - self.r_bins[i]
                dvol = 2*np.arccos(-1)*self.r_bins[i]*dr
            if norm_type == '3d':
                dr = self.r_bins[i+1] - self.r_bins[i]
                dvol = 4*np.arccos(-1)*(self.r_bins[i]**2)*dr

            self.histogram[i]/=dvol*norm_factor

    
    def getHistogram(self):
        out_hist = []
        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            out_hist.append( [ r, self.histogram[i] ] )
        return out_hist

    def print_to(self, stream):
        headline='# '+str(self.r_bin_nb)+' '+str(self.r_max)+'\n'
        stream.write(headline)

        for i in range(self.r_bin_nb):
            r=self.r_bins[i]
            out_string=str(r)+" "+str(self.histogram[i])+'\n'
            stream.write(out_string)


    
