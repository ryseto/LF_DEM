#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import math
import pycart2sph
import numpy as np
import string

class SphericalCoordinateHistogram:

    def __init__(self, params):
        [r_bn, theta_bn, phi_bn, min_r, max_r, min_theta, max_theta] = params
        
        self.r_bin_nb=r_bn
        self.theta_bin_nb=theta_bn
        self.phi_bin_nb=phi_bn

        self.r_max=max_r
        self.r_min=min_r
        self.theta_max=max_theta
        self.theta_min=min_theta
        self.phi_min=0 # for future extension

        self.pi=np.arccos(-1)
        self.r_bsize=(self.r_max-self.r_min)/self.r_bin_nb
        self.theta_bsize=(self.theta_max-self.theta_min)/self.theta_bin_nb
        self.phi_bsize=self.pi/self.phi_bin_nb

        self.histogram = np.zeros((r_bn, theta_bn, phi_bn))



    def update(self, coord, value, coord_sys='sph', act='add'):
    
        if coord_sys == 'cart':
            (r, theta, phi)=pycart2sph.cart2sph(coord)
        else:
            coordarray = np.array(coord)
            r = coordarray[:,0]
            theta = coordarray[:,1]
            phi = coordarray[:,2]

        r_bin = np.array((r-self.r_min)/self.r_bsize, dtype=int)
        theta_bin = np.array((theta-self.theta_min)/self.theta_bsize, dtype=int)
        phi_bin = np.array((phi-self.phi_min)/self.phi_bsize, dtype=int)
        

        inrange_bin = np.logical_and(r_bin<self.r_bin_nb, r_bin>=0)
        inrange_bin = np.logical_and(inrange_bin, np.logical_and(theta_bin<self.theta_bin_nb, theta_bin>=0))
        inrange_bin = np.logical_and(inrange_bin, np.logical_and(phi_bin<self.phi_bin_nb, phi_bin>=0))

        r_bin = r_bin[ inrange_bin ]
        theta_bin = theta_bin[ inrange_bin ]
        phi_bin = phi_bin[ inrange_bin ]
        values = np.array(value)[ inrange_bin ]

        
        if act == 'add':
            # 'add' is a bit tricky, 
            # because treats a[indices] += x in a somewhat unexpected way:
            # it ignores duplicates in the indices
            # the solution here is from:
            # http://stackoverflow.com/questions/18273634/assigning-identical-array-indices-at-once-in-python-numpy
            # the idea is to use bincount to aggregate the duplicate indices first, and then feed the array with an array of unique indices. 
            # but, as things cannot be simple, bincount can only aggregate 1d indices, 
            # so there's a pre-step to convert 3d indices in 1d, with ravel_multi_index
            flat_indices = np.ravel_multi_index((r_bin, theta_bin, phi_bin), dims=self.histogram.shape)
            unique_indices = np.unique(flat_indices)
            aggregated_values = (np.bincount(flat_indices, weights = values))[unique_indices]
            self.histogram[np.unravel_index(unique_indices, dims=self.histogram.shape)] += aggregated_values
        elif act == 'replace':
            self.histogram[r_bin, theta_bin, phi_bin] = values
        else:
            sys.stderr.write("  SphericalCoordinateHistogram :: update(pos, value, coord='sph',act='add') :\n       Unknown act %s".str(act))
            sys.exit(1)


        
    def normalize(self, norm_factor=1): # for g(r), norm_factor = V*nb_of_snapshots
        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            rlo=self.r_bsize*(i+0.001)+self.r_min
            for j in range(self.theta_bin_nb):
                theta=self.theta_bsize*(j+0.5)+self.theta_min
                for k in range(self.phi_bin_nb):
                    phi=self.phi_bsize*(k+0.5)+self.phi_min
                    dvol=rlo*rlo*math.sin(theta)*self.theta_bsize*self.phi_bsize*self.r_bsize
                    self.histogram[i,j,k]/=dvol*norm_factor

    # def bin_coord(self,key):
    #     [i,j,k]=string.split(key)        
    #     r=self.r_bsize*(float(i)+0.5)+self.r_min
    #     theta=self.theta_bsize*(float(j)+0.5)
    #     phi=self.phi_bsize*(float(k)+0.5)
    #     return (r, theta,phi)
    


    def print_to(self, stream): # to be vectorized...
        headline=str(self.r_bin_nb)+' '+str(self.theta_bin_nb)+' '+str(self.phi_bin_nb)+' '+str(self.r_max)+'\n'
        stream.write(headline)

        for i in range(self.r_bin_nb):
            r=self.r_bsize*(i+0.5)+self.r_min
            for j in range(self.theta_bin_nb):
                theta=self.theta_bsize*(j+0.5)+self.theta_min
                for k in range(self.phi_bin_nb):
                    phi=self.phi_bsize*(k+0.5)+self.phi_min
#                out_string=str(r)+" "+str(theta)+" "+str(phi)+" "+" ".join(map(str, self.histogram[point]))+'\n'
                    out_string=str(r)+" "+str(theta)+" "+str(phi)+" "+str(self.histogram[i,j,k])+'\n'
                    stream.write(out_string)


    
