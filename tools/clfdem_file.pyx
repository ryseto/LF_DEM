# distutils: language = c++
# distutils: sources = lfdem_files.cpp

from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp cimport bool

import numpy as np


cdef extern from "lfdem_files.cpp":
    cdef struct Frame:
        map[string, double] meta_data
        vector[vector[double]] data

    cdef cppclass lf_file:
        lf_file(string) except +
        map[string, string] get_meta_data()
        map[string, pair [int, int]] get_column_def()

    cdef cppclass lf_data_file(lf_file):
        lf_data_file(string) except +
        vector[vector[double]] get_data()

    cdef cppclass lf_snapshot_file(lf_file):
        lf_snapshot_file(string) except +
        Frame get_frame(unsigned int)
        Frame next_frame()


cdef class generic_file:
    cdef lf_file *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, fname):
        if type(self) is generic_file:
            self.thisptr = new lf_file(fname.encode('utf8'))

    def __dealloc__(self):
        if type(self) is generic_file:
            del self.thisptr

    def meta_data(self):
        return self.thisptr.get_meta_data()

    def column_def(self):
        return self.thisptr.get_column_def()

cdef class data_file(generic_file):
    cdef lf_data_file *derivedthisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, fname):
        self.derivedthisptr = self.thisptr = new lf_data_file(fname.encode('utf8'))

    def __dealloc__(self):
        del self.derivedthisptr

    def data(self):
      return np.array(self.derivedthisptr.get_data(), dtype=np.float)

cdef class snapshot_file(generic_file):
    cdef lf_snapshot_file *derivedthisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, fname):
        self.derivedthisptr = self.thisptr = new lf_snapshot_file(fname.encode('utf8'))

    def __dealloc__(self):
        del self.derivedthisptr

    def next_frame(self):
        frame = self.derivedthisptr.next_frame()
        return frame.meta_data, np.array(frame.data, dtype=np.float)
