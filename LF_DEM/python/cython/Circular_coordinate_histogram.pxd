#coding=utf-8

cimport Histograms
from Histograms cimport Histograms

cdef class CircularCoordinateHistogram(Histograms):
    cdef:
        int r_bin_nb
        int theta_bin_nb
        double r_max
        double r_min
        double r_bsize
        double theta_bsize
        double pi

    cpdef genkey_nolist(self, double r, double theta)
    cpdef update_nolist(self, double a, double b, value, str coord=*, str act=*)

    
