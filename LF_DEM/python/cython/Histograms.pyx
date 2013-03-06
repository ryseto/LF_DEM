#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import math
import numpy
cimport numpy
import string
import copy


cdef class Histograms:

    def __init__(self):
        self.histogram=dict()
        self.events=dict()

    def __iter__(self):
        return self.histogram.iterkeys()

    def copy(self, histo):
        self.histogram=histo.histogram_copy()
        self.events=histo.events_copy()

    def deepcopy(self, histo):
        self.histogram=histo.histogram_deepcopy()
        self.events=histo.events_deepcopy()

    cdef kadd(self, key, value):
        if key is not None:
            value_array=numpy.array(value)
            if self.histogram[key] is not None:
                self.histogram[key]+=value_array
            else:
                self.histogram[key]=value_array
            self.events[key]+=1
            

    cdef kappend(self, key, value):
        if key is not None:
            value_array=numpy.array([value])
            if self.histogram[key] is not None:
                new_array=numpy.append(self.histogram[key], value_array, axis=0)
                self.kreplace(key, new_array)
                self.events[key]+=1
            else:
                self.histogram[key]=value_array

    cdef kreplace(self, key, value):
        if key is not None:
            value_array=numpy.array(value)
            self.histogram[key]=value_array

    def normalize(self):
        for point in self.histogram:
            if self.histogram[point] is not None:
                self.histogram[point]/=float(self.events[point])

    def print_to(self, stream):
        for point in self.histogram:
            out_string=str(point)+" "+" ".join(map(str, self.histogram[point]))+'\n'

    def histogram_copy(self):
        return self.histogram.copy()

    def events_copy(self):
        return self.events.copy()

    def histogram_deepcopy(self):
        return copy.deepcopy(self.histogram)

    def events_deepcopy(self):
        return copy.deepcopy(self.events)

    def kvalue(self, key):
        return self.histogram[key]

    
