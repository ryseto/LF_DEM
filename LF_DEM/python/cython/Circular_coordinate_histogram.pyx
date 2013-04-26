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

cdef class CircularCoordinateHistogram(Histograms):
    # cdef:
    #     int r_bin_nb
    #     int theta_bin_nb
    #     double r_max
    #     double r_min
    #     double r_bsize
    #     double theta_bsize
    #     double pi

    def __init__(self, int r_bn, int theta_bn, double min_r, double max_r):
        Histograms.__init__(self)
        
        self.r_bin_nb=r_bn
        self.theta_bin_nb=theta_bn

        self.r_max=max_r
        self.r_min=min_r

        self.pi=3.141592653589793238462643
        self.r_bsize=(max_r-self.r_min)/self.r_bin_nb
        self.theta_bsize=self.pi/self.theta_bin_nb # symmetry: only theta in [0:pi]

        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            for j in range(self.theta_bin_nb):
                theta=self.theta_bsize*(j+0.5)
                key=self.genkey_nolist(r, theta)
                self.histogram[key]=None
                self.events[key]=0


    def __iter__(self):
        return self.histogram.iterkeys()

    def copy(self, histo):
        Histograms.copy(histo)
        (self.r_bin_nb, self.theta_bin_nb, self.r_max, self.r_bsize, self.theta_bsize)=histo.params()

    def deepcopy(self, histo):
        Histograms.deepcopy(histo)
        (self.r_bin_nb, self.theta_bin_nb, self.r_max, self.r_bsize, self.theta_bsize)=histo.params()

    cpdef genkey_nolist(self, double r, double theta):
        cdef int r_bin_label
        cdef int theta_bin_label

        #take advantage of the symmetry of the flow
        if theta>self.pi:
            theta = theta - self.pi

        r_bin_label=int((r-self.r_min)/self.r_bsize)
        theta_bin_label=int(theta/self.theta_bsize)

        if theta_bin_label==self.theta_bin_nb:
            theta_bin_label-=self.theta_bin_nb

        if r_bin_label>=self.r_bin_nb:
            return None


        key=str(r_bin_label)+' '+str(theta_bin_label)
        return key


    def genkey(self, sph_pos):
        return self.genkey_nolist(sph_pos[0], sph_pos[1])


    cpdef update_nolist(self, double a, double b, value, str coord='circ', str act='add'):

        if coord == 'cart':
            (a, b)=cart2sph.cart2circ_nolist(a, b)

        cdef str key
        key=self.genkey_nolist(a, b)

        if act == 'append':
            self.kappend(key, value)
        elif act == 'add':
            self.kadd(key, value)
        elif act == 'replace':
            self.kreplace(key, value)
        else:
            sys.stderr.write("  CircularCoordinateHistogram :: update(pos, value, coord='sph',act='add') :\n       Unknown act %s".str(act))
            sys.exit(1)


    def update(self, pos, value, coordinates='circ',action='add'):
        self.update_nolist(pos[0], pos[1], value, coord=coordinates,act=action)
        
    def normalize_a_la_gofr(self, int N, double rho, int call_nb):

        for point in self.histogram:
            if self.histogram[point] is None:
                self.histogram[point] = 0.

        norm_factor=N*rho*call_nb
        
        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            rlo=self.r_bsize*i+self.r_min
            for j in range(self.theta_bin_nb):
                theta=self.theta_bsize*(j+0.5)
                key=self.genkey([r, theta])
                dvol=rlo*self.theta_bsize*self.r_bsize
                self.histogram[key]/=dvol*norm_factor

    def bin_coord(self,key):
        [i,j]=string.split(key)
        r=self.r_bsize*(float(i)+0.5)+self.r_min
        theta=self.theta_bsize*(float(j)+0.5)
        return (r, theta)
    
    def print_to(self, stream):
        headline=str(self.r_bin_nb)+' '+str(self.theta_bin_nb)+' '+str(self.r_max)+'\n'
        stream.write(headline)

        for point in self.histogram:
            if self.histogram[point] is not None:
                (r, theta)=self.bin_coord(point)
#                out_string=str(r)+" "+str(theta)+" "+str(phi)+" "+" ".join(map(str, self.histogram[point]))+'\n'
                out_string=str(r)+" "+str(theta)+" "+str(self.histogram[point])+'\n'
                stream.write(out_string)

    def params(self):
        return (self.r_bin_nb, self.theta_bin_nb, self.r_max, self.r_bsize, self.theta_bsize)


    
