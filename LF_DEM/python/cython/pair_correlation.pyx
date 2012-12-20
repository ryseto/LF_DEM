#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import LF_DEM_posfile_reading
import Spherical_coordinate_histogram
import cart2sph

cdef class PairCorrelation:

    cdef:
        field
    def __init__(self, r_bin_nb,theta_bin_nb,phi_bin_nb, r_max):
        self.field=Spherical_coordinate_histogram.SphericalCoordinateHistogram(r_bin_nb,theta_bin_nb,phi_bin_nb, r_max)


    cpdef update_field(self,pos_stream):
        cdef int i
        cdef int j
        cdef double deltax, deltay, deltaz, deltatot

        for i in pos_stream.range():
            for j in pos_stream.range(i+1):
                #                pos_stream.cpos_velocity_diff(i,j, &deltax , &deltay, &deltaz, &deltatot, &velx, &vely, &velz, &veltot)

                dr=pos_stream.pos_diff(i,j)
                print "dr " , dr
                deltax=dr[0]
                deltay=dr[1]
                deltaz=dr[2]
                deltatot=dr[3]
                
                (a, b, c)=cart2sph.cart2sph_nolist(deltax, deltay, deltaz)
                print "dr_sph " , [a, b, c]
                key=self.field.genkey_nolist(a, b, c)
                
                self.field.update_nolist(a, b, c, 1., coord='sph', act='add')

    def normalize(self):
        self.field.normalize()

    def print_to(self, stream):
        self.field.print_to(stream)
