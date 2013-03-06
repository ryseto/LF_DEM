#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import LF_DEM_posfile_reading
cimport Spherical_coordinate_histogram
from Spherical_coordinate_histogram cimport SphericalCoordinateHistogram
cimport Circular_coordinate_histogram

import cart2sph

cdef class PairCorrelation:

    cdef:
        field
        int update_nb
        int dimension
    def __init__(self, dim, params ):

        if dim==3:
            self.field=SphericalCoordinateHistogram(params)
        if dim==2:
            [r_bn, theta_bn, min_r, max_r] = params
            self.field=Circular_coordinate_histogram.CircularCoordinateHistogram(r_bn, theta_bn, min_r, max_r)
        self.dimension = dim


    cdef update_field_2d(self,pos_stream):

        cdef int i
        cdef int j
        cdef double deltax, deltay, deltaz, deltatot

        
        self.update_nb += 1

        for i in pos_stream.range():

            for j in pos_stream.range(i+1):
                dr=pos_stream.pos_diff(i,j)

                deltax=dr[0]
                deltaz=dr[2]
                deltatot=dr[3]

                # we go to circular coordinates
                (a, b)=cart2sph.cart2circ_nolist(deltax, deltaz)

                key=self.field.genkey_nolist(a, b)
                self.field.update_nolist(a, b, 1., coord='circ', act='add')



    cdef update_field_3d(self,pos_stream):
        cdef int i
        cdef int j
        cdef double deltax, deltay, deltaz, deltatot

        
        self.update_nb += 1

        for i in pos_stream.range():

            for j in pos_stream.range(i+1):
                dr=pos_stream.pos_diff(i,j)

                deltax=dr[0]
                deltay=dr[1]
                deltaz=dr[2]
                deltatot=dr[3]

                # we go to spherical coordinates
                # note that z and y directions are exchanged: it is not a bug
                # theta=pi/2. is the shear plane
                # theta=0 is the vorticity direction (y)
                (a, b, c)=cart2sph.cart2sph_nolist(deltax, deltaz, deltay)

                key=self.field.genkey_nolist(a, b, c)
                self.field.update_nolist(a, b, c, 1., coord='sph', act='add')





    cpdef update_field(self,pos_stream):

        if self.dimension == 2:
            self.update_field_2d(pos_stream)
        if self.dimension == 3:
            self.update_field_3d(pos_stream)




    def normalize(self, N, rho):
        self.field.normalize_a_la_gofr(N, rho, self.update_nb)

    def print_to(self, stream):
        self.field.print_to(stream)
