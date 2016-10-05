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
        void rewind()

    cdef cppclass lf_data_file(lf_file):
        lf_data_file(string) except +
        vector[vector[double]] get_data()

    cdef cppclass lf_snapshot_file(lf_file):
        lf_snapshot_file(string) except +
        Frame get_frame(size_t)
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
        coldef_dict = dict(self.thisptr.get_column_def())
        coldef_dict = dict(zip([a.decode() for a in coldef_dict], coldef_dict.values()))
        for k in coldef_dict:
            if coldef_dict[k][1]>coldef_dict[k][0]:
                colslice = slice(coldef_dict[k][0]-1, coldef_dict[k][1])
                coldef_dict[k] = colslice
            else:
                coldef_dict[k] = coldef_dict[k][0]
        return coldef_dict

    def rewind(self):
        return self.thisptr.rewind()

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

    def __iter__(self):
        self.derivedthisptr.rewind()
        return self

    def __next__(self):
        frame = self.derivedthisptr.next_frame()
        if len(frame.meta_data):
            meta_dict = dict(frame.meta_data)
            meta_dict = dict(zip([a.decode() for a in meta_dict], meta_dict.values()))
            return frame.meta_data, np.array(frame.data, dtype=np.float)
        else:
            raise StopIteration

    def __getitem__(self, index):
        frame = self.derivedthisptr.get_frame(index)
        if len(frame.meta_data):
            return frame.meta_data, np.array(frame.data, dtype=np.float)
        else:
            raise IndexError
