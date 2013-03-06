#coding=utf-8

cimport Histograms
from Histograms cimport Histograms

cdef class LinearHistogram(Histograms):
    cdef:
        int r_bin_nb
        double r_max
        double r_min
        double r_bsize

    cpdef genkey_nolist(self, double r)
    cpdef update_nolist(self, double a, value, str act=*)

    
