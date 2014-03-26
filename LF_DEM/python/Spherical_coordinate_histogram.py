#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import math
import cart2sph
import numpy as np
import string

class SphericalCoordinateHistogram:

    def __init__(self, params):
        self.pi=np.arccos(-1)
        if params[0] == 'linear':
            [ r_bn, theta_bn, phi_bn, min_r, max_r, min_theta, max_theta ] = params[1:]
            self.r_bins = np.linspace(self.min_r, self.max_r, self.r_bn)
            self.theta_bins = np.linspace(self.min_theta, self.max_theta, self.theta_bn)
            self.phi_bins = np.linspace(0, self.pi, self.phi_bin_nb)
        elif params[0] == 'custom':
            [ self.r_bins, self.theta_bins, self.phi_bins ] = [ np.array(params[1]), np.array(params[2]), np.array(params[3]) ]
            

        self.r_bin_nb = self.r_bins.shape[0]-1
        self.r_max = self.r_bins.max()
        self.r_min = self.r_bins.min()
        self.r_bsize=(self.r_max-self.r_min)/self.r_bin_nb

        self.theta_bin_nb = self.theta_bins.shape[0]-1
        self.theta_max = self.theta_bins.max()
        self.theta_min = self.theta_bins.min()
        self.theta_bsize=(self.theta_max-self.theta_min)/self.theta_bin_nb

        self.phi_bin_nb = self.phi_bins.shape[0]-1
        self.phi_max = self.phi_bins.max()
        self.phi_min = self.phi_bins.min()
        self.phi_bsize=self.pi/self.phi_bin_nb

        self.histogram = np.zeros((self.r_bin_nb, self.theta_bin_nb, self.phi_bin_nb))


    def update(self, coord, value, coord_sys='sph', act='add'):
        
        if coord_sys == 'cart':
            (r, theta, phi)=pycart2sph.cart2sph(coord)
        else:
            coordarray = np.array(coord)
            r = coordarray[:,0]
            theta = coordarray[:,1]
            phi = coordarray[:,2]


        inrange = np.logical_and(r>self.r_min, r<self.r_max)
        inrange = np.logical_and(inrange, np.logical_and(theta>self.theta_min, theta<self.theta_max))
        inrange = np.logical_and(inrange, np.logical_and(phi>self.phi_min, phi<self.phi_max))
        inrange = np.nonzero(inrange)

        if len(inrange[0]) == 0:
            return

        r = r [ inrange ]
        theta = theta [ inrange ]
        phi = phi [ inrange ]
        values = np.array(value)[ inrange ]
        
#        print r, theta, phi, values
        r_bin = np.digitize(r,self.r_bins)-1
        theta_bin = np.digitize(theta,self.theta_bins)-1
        phi_bin = np.digitize(phi,self.phi_bins)-1

        
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
            r = self.r_bins[i]
            dr = self.r_bins[i+1] - self.r_bins[i]
            for j in range(self.theta_bin_nb):
                theta = 0.5*(self.theta_bins[j+1] + self.theta_bins[j])
                dtheta = self.theta_bins[j+1] - self.theta_bins[j]
                for k in range(self.phi_bin_nb):
                    dphi = self.phi_bins[k+1] - self.phi_bins[k] 

                    dvol= (r**2)*dr*np.sin(theta)*dtheta*dphi
                    self.histogram[i,j,k]/=dvol*norm_factor

    def getHistogram(self):
#        out_hist = np.zeros((r_bin_nb*theta_bin_nb*phi_bin_nb))
        out_hist =[]
        for i in range(self.r_bin_nb):
            r=self.r_bins[i]
            for j in range(self.theta_bin_nb):
                theta=self.theta_bsize*(j+0.5)+self.theta_min
                for k in range(self.phi_bin_nb):
                    phi=self.phi_bsize*(k+0.5)+self.phi_min
                    out_hist.append( [ r, theta, phi, self.histogram[i,j,k] ] )
        return np.array(out_hist)

    def print_to(self, stream): # to be vectorized...
        headline='# '+str(self.r_bin_nb)+' '+str(self.theta_bin_nb)+' '+str(self.phi_bin_nb)+' '+str(self.r_max)+'\n'
        stream.write(headline)

        for i in range(self.r_bin_nb):
            r=self.r_bins[i]
            for j in range(self.theta_bin_nb):
                theta=self.theta_bsize*(j+0.5)+self.theta_min
                for k in range(self.phi_bin_nb):
                    phi=self.phi_bsize*(k+0.5)+self.phi_min
#                out_string=str(r)+" "+str(theta)+" "+str(phi)+" "+" ".join(map(str, self.histogram[point]))+'\n'
                    out_string=str(r)+" "+str(theta)+" "+str(phi)+" "+str(self.histogram[i,j,k])+'\n'
                    stream.write(out_string)


    
