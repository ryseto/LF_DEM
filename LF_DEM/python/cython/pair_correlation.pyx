#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

from __future__ import division
import sys
import math
import pyLF_DEM_posfile_reading
import Spherical_coordinate_histogram
from Spherical_coordinate_histogram import SphericalCoordinateHistogram
import Circular_coordinate_histogram
from Circular_coordinate_histogram import CircularCoordinateHistogram
import Linear_histogram
from Linear_histogram import LinearHistogram

import cart2sph
import numpy as np

cdef class PairCorrelation:

    cdef:
        field
        int update_nb
        int mode
        restriction_couples
        max_out_of_plane
        rescale_func

    def __init__(self, _mode, params, restrictions=None, rescale=None):

        # restrictions have the form:
        # ( function, value1, value2)
        # pairs of particles (i,j) taken into account in the partition function
        # verify ( function(i)==value1 and function(j)==value2 ) or ( function(j)==value1 and function(i)==value2 )
        # several restrictions can be given, in which case all of them have to be satisfied
        
        self.restriction_couples = list(restrictions) if restrictions is not None else []
        self.rescale_func = rescale if rescale is not None else lambda x:0.5

        if _mode == "r":
            self.mode=0
            [r_bn, min_r, max_r] = params
            self.field=LinearHistogram(r_bn, min_r, max_r)

        if _mode == "c":
            self.mode=1
            self.max_out_of_plane = 0.05
            if len(params) == 5:
                [r_bn, theta_bn, min_r, max_r, self.max_out_of_plane] = list(params) # for 3d systems
            elif len(params) == 4:
                [r_bn, theta_bn, min_r, max_r ] = list(params)
            else:
                sys.stderr.write('Circular pair correlation init : invalid list of parameters')

            self.field=CircularCoordinateHistogram(r_bn, theta_bn, min_r, max_r)

        if _mode == "s":
            self.mode=2
            self.field=SphericalCoordinateHistogram(params)            


    cdef update_field_1d(self, double deltatot):
    
        self.field.update_nolist(deltatot, 1., act='add')



    cdef update_field_2d(self,double deltax, double deltay, double deltaz, double deltatot):

        if deltatot>0 and math.fabs(deltay)/deltatot < self.max_out_of_plane:
                
            # we go to circular coordinates
            (a, b)=cart2sph.cart2circ_nolist(deltax, deltaz)

            self.field.update_nolist(a, b, 1., coord='circ', act='add')



    cdef update_field_3d(self, double deltax, double deltay, double deltaz):

        # we go to spherical coordinates
        # note that z and y directions are exchanged: it is not a bug
        # theta=pi/2. is the shear plane
        # theta=0 is the vorticity direction (y)
        (a, b, c)=cart2sph.cart2sph_nolist(deltax, deltaz, deltay)

        self.field.update_nolist(a, b, c, 1., coord='sph', act='add')
    

    cdef passRestrictions(self, i, j):
    
        for r in self.restriction_couples:
            if not ( ( r[0](i) == r[1] and r[0](j) == r[2] ) or ( r[0](j) == r[1] and r[0](i) == r[2] ) ):
                return False

        return True

    cpdef update_field(self,pos_stream):

        cdef int i
        cdef int j
        cdef double deltax, deltay, deltaz, deltatot

        self.update_nb += 1
        
        pos_stream.computePairSeparations()
        for i in pos_stream.range():
            for j in pos_stream.range(i+1):
                if self.passRestrictions(i, j):
                    rescale_factor = 1/(self.rescale_func(i) + self.rescale_func(j))
                    deltax = pos_stream.pair_sep[i,j-i][0]*rescale_factor
                    deltay = pos_stream.pair_sep[i,j-i][1]*rescale_factor
                    deltaz = pos_stream.pair_sep[i,j-i][2]*rescale_factor
                    deltatot = pos_stream.pair_dist[i,j-i]*rescale_factor
                    
                    if self.mode == 0:
                        self.update_field_1d(deltatot)
                    if self.mode == 1:
                        self.update_field_2d(deltax, deltay, deltaz, deltatot)
                    if self.mode == 2:
                        self.update_field_3d(deltax, deltay, deltaz)
                                


    def normalize(self, Nrho):
        self.field.normalize_a_la_gofr(Nrho, self.update_nb)

    def print_to(self, stream):
        self.field.print_to(stream)
    
    def getField(self):
        g_of_r = []
        for point in self.field.histogram_copy().iteritems():
            g_of_r.append( [ self.field.bin_coord(point[0]), point[1] ] )
        return  g_of_r

    def getUpdateNb(self):
        return self.update_nb
