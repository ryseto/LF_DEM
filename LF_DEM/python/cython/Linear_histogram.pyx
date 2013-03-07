#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import math
import cart2sph
import numpy
cimport numpy
import string
import copy
cimport Histograms
from Histograms cimport Histograms

cdef class LinearHistogram(Histograms):

    def __init__(self, int r_bn, double min_r, double max_r):
        Histograms.__init__(self)
        
        self.r_bin_nb=r_bn

        self.r_max=max_r
        self.r_min=min_r

        self.r_bsize=(max_r-self.r_min)/self.r_bin_nb

        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            key=self.genkey_nolist(r)
            self.histogram[key]=None
            self.events[key]=0


    def __iter__(self):
        return self.histogram.iterkeys()

    def copy(self, histo):
        Histograms.copy(histo)
        (self.r_bin_nb, self.r_max, self.r_bsize, self.theta_bsize)=histo.params()

    def deepcopy(self, histo):
        Histograms.deepcopy(histo)
        (self.r_bin_nb, self.r_max, self.r_bsize, self.theta_bsize)=histo.params()

    cpdef genkey_nolist(self, double r):
        cdef int r_bin_label

        r_bin_label=int((r-self.r_min)/self.r_bsize)

        if r_bin_label>=self.r_bin_nb:
            return None

        key=str(r_bin_label)
        return key


    def genkey(self, coord):
        return self.genkey_nolist(coord)


    cpdef update_nolist(self, double a, value, str act='add'):

        cdef str key
        key=self.genkey_nolist(a)

        if act == 'append':
            self.kappend(key, value)
        elif act == 'add':
            self.kadd(key, value)
        elif act == 'replace':
            self.kreplace(key, value)
        else:
            sys.stderr.write("  SphericalCoordinateHistogram :: update(pos, value, coord='sph',act='add') :\n       Unknown act %s".str(act))
            sys.exit(1)


    def update(self, coord, value, action='add'):
        self.update_nolist(coord, value, act=action)
        
    def normalize_a_la_gofr(self, int N, double rho, int call_nb):

        for point in self.histogram:
            if self.histogram[point] is None:
                self.histogram[point] = 0.

        norm_factor=N*rho*call_nb
        
        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            rlo=self.r_bsize*i+self.r_min
            key=self.genkey(r)
            dvol=rlo*self.r_bsize
            self.histogram[key]/=dvol*norm_factor

    def bin_coord(self,key):
        i=int(key)
        r=self.r_bsize*(float(i)+0.5)+self.r_min
        return r
    
    def print_to(self, stream):
        headline=str(self.r_bin_nb)+' '+str(self.r_max)+'\n'
        stream.write(headline)

        for point in self.histogram:
            if self.histogram[point] is not None:
                r=self.bin_coord(point)
                out_string=str(r)+" "+str(self.histogram[point])+'\n'
                stream.write(out_string)

    def params(self):
        return (self.r_bin_nb, self.r_max, self.r_bsize, self.theta_bsize)


    
