# distutils: language = c++
# distutils: sources = lfdem_files.cpp

from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp cimport bool

import numpy as np
import struct as cstruct


cdef extern from "lfdem_files.cpp":
    cdef struct Frame:
        map[string, double] meta_data
        vector[vector[double]] data

    cdef cppclass lf_file:
        lf_file(string) except +
#        lf_file(string, vector[int]) except +
        lf_file(string, vector[string]) except +

        map[string, string] get_meta_data()
        map[string, pair [int, int]] get_column_def()
        void rewind()

    cdef cppclass lf_data_file(lf_file):
        lf_data_file(string) except +
#        lf_data_file(string, vector[int]) except +
        lf_data_file(string, vector[string]) except +

        vector[vector[double]] get_data()

    cdef cppclass lf_snapshot_file(lf_file):
        lf_snapshot_file(string) except +
        Frame get_frame(size_t) except +
        Frame next_frame() except +


cdef class generic_file:
    cdef lf_file *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, fname, fields_to_read=None):
        if type(self) is generic_file:
            if fields_to_read is None:
                self.thisptr = new lf_file(fname.encode('utf8'))
            else:
                self.thisptr = new lf_file(fname.encode('utf8'),
                                           fields_to_read)


    def __dealloc__(self):
        if type(self) is generic_file:
            del self.thisptr

    def meta_data(self):
        meta_dict = dict(self.thisptr.get_meta_data())
        meta_dict = dict(zip([a.decode() for a in meta_dict],
                             [a.decode() for a in meta_dict.values()]))
        return meta_dict

    def column_def(self):
        coldef_dict = dict(self.thisptr.get_column_def())
        coldef_dict = dict(zip([a.decode() for a in coldef_dict],
                               coldef_dict.values()))
        for k in coldef_dict:
            if coldef_dict[k][1] > coldef_dict[k][0]+1:
                colslice = slice(coldef_dict[k][0], coldef_dict[k][1])
                coldef_dict[k] = colslice
            else:
                coldef_dict[k] = coldef_dict[k][0]
        return coldef_dict

    def rewind(self):
        return self.thisptr.rewind()

cdef class data_file(generic_file):
    r"""
    LF_DEM data file object, i.e. files with one line per record (e.g. data_ or st_ files)

    Arguments:
        fname: file name
    Optional argument:
        fields_to_read: iterable of column names, e.g. ["viscosity", "shear strain"] (default is all fields)
    """

    cdef lf_data_file *derivedthisptr      # hold a C++ instance which we're wrapping

    def __cinit__(self, fname, fields_to_read=None):
        if fields_to_read is None:
            self.derivedthisptr = self.thisptr = new lf_data_file(fname.encode('utf8'))
        else:
            if isinstance(fields_to_read[0], str):
                fields_to_read = [c.encode('utf8') for c in fields_to_read]
                self.derivedthisptr = self.thisptr = new lf_data_file(fname.encode('utf8'),
                                                                      fields_to_read)

    def __dealloc__(self):
        del self.derivedthisptr

    def data(self):
        """Returns a numpy array of the data read from file"""
        return np.array(self.derivedthisptr.get_data(), dtype=np.float)

cdef class snapshot_file(generic_file):
    r"""
    LF_DEM snapshot file object, i.e. files with multiple lines per record (e.g. par_ or int_ files)
    This object is iterable, which means you can iterate over records.

    Arguments:
        fname: file name
    """

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
            return meta_dict, np.array(frame.data, dtype=np.float)
        else:
            raise StopIteration

    def __getitem__(self, index):
        frame = self.derivedthisptr.get_frame(index)
        if len(frame.meta_data):
            return frame.meta_data, np.array(frame.data, dtype=np.float)
        else:
            raise IndexError

def column_def(fname):
    f = generic_file(fname)
    return f.column_def()


def read_conf_file(fname):
    """
    Purpose:
        Read any LF_DEM conf file (usually D?N?VF?_*.dat files)

    Parameters:
        fname: the filename, or anything that can be taken as a first argument
               to np.genfromtxt

    Returning values:
        positions, radii, metadata
    """
    openedfile = True
    try:
        in_file = open(fname, "r")
    except TypeError:
        in_file = fname
        openedfile = False

    line1 = in_file.readline()
    line2 = in_file.readline()
    meta_fields = line1.split()[1:]
    meta_values = line2.split()[1:]
    meta_data = dict(zip(meta_fields, meta_values))

    dat = np.genfromtxt(fname, usecols=np.arange(4))
    pos = dat[:, :3].astype(np.float)
    rad = dat[:, 3].astype(np.float)
    if openedfile:
        in_file.close()
    return pos, rad, meta_data


def popValue(t, stream):
    size = cstruct.calcsize(t)
    buf = stream.read(size)
    return cstruct.unpack(t, buf)


def read_binary_conf_file(fname):

    stream = open(fname, mode='rb')
    meta_data = {}
    config = {}

    # determine the format
    i = popValue("i", stream)[0]
    if i == -1:
        meta_data["format"], meta_data["np"] = popValue("2i", stream)
        if meta_data["format"] > 3:
            print(" unknown LF_DEM binary format : ", meta_data["format"])
            exit(1)
        if meta_data["format"] == 3:
            meta_data["np_fixed"] = popValue("i", stream)[0]
    else:
        meta_data["np"] = i
        meta_data["format"] = 2
    meta_data["vf"] = popValue("d", stream)[0]
    meta_data["lx"], meta_data["ly"], meta_data["lz"] = popValue("3d", stream)
    meta_data["disp_x"], meta_data["disp_y"] = popValue("2d", stream)
    config['metadata'] = meta_data
    config['positions'] = np.empty((meta_data["np"], 4))
    for i in range(meta_data["np"]):
        config['positions'][i] = np.array(popValue("4d", stream))

    if meta_data["format"] == 3:
        config['fixed_velocities'] = np.empty((meta_data["np_fixed"], 3))
        for i in range(meta_data["np_fixed"]):
            config['fixed_velocities'][i] = np.array(popValue("3d", stream))

    nc = popValue("I", stream)[0]
    config['contacts'] = np.empty((nc, 8))
    for i in range(nc):
        config['contacts'][i] = np.array(popValue("2I6d", stream))

    stream.close()
    return config


def pushValue(stream, fmt, *values):
    stream.write(cstruct.pack(fmt, *values))


def write_binary_conf_file(fname,
                           config):
    mandatory_fields = ['metadata', 'positions', 'contacts']
    for field in mandatory_fields:
        if field not in config:
            raise RuntimeError("Config must contain a "+field+" field.")

    with open(fname, "wb") as f:
        meta_data = config["metadata"]
        pushValue(f, "3i", -1, meta_data["format"], meta_data["np"])
        if meta_data["format"] == 3:
            pushValue(f, "i", meta_data["np_fixed"])
        pushValue(f, "6d",
                  meta_data["vf"],
                  meta_data["lx"],
                  meta_data["ly"],
                  meta_data["lz"],
                  meta_data["disp_x"],
                  meta_data["disp_y"])

        pos_size = str(config["positions"].size)
        pushValue(f, pos_size+"d", *np.ravel(config["positions"]))

        if meta_data["format"] == 3:
            vel_size = str(config["fixed_velocities"].size)
            pushValue(f, vel_size+"d", *np.ravel(config["fixed_velocities"]))

        pushValue(f, "I", config['contacts'].shape[0])
        for contact in config["contacts"]:
            pushValue(f, "2I",*contact[:2].astype(np.int))
            pushValue(f, "6d",*contact[2:])
