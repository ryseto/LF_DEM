#coding=utf-8

cimport Histograms
from Histograms cimport Histograms

cdef class SphericalCoordinateHistogram(Histograms):
    cdef:
        int r_bin_nb
        int theta_bin_nb
        int phi_bin_nb
        double r_max
        double r_min
        double r_bsize
        double theta_bsize
        double phi_bsize
        double pi

    cpdef genkey_nolist(self, double r, double theta, double phi)
    cpdef update_nolist(self, double a, double b, double c, value, str coord=*, str act=*)
#    cpdef normalize_a_la_gofr(self, int N, double rho, int call_nb)
