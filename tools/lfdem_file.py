# lfdem_file.py, part of LF_DEM
#
#  Romain Mari, 2015
"""lfdem_file is a module providing functions to read output files from lf_dem"""

import pandas as pd
import numpy as np



def read_snapshot_file(in_file, field_nb=None):
    """
    Purpose:
        Read any LF_DEM file that has a "snapshot" structure, i.e. made of the following:
            a header (x lines starting with #)
            an empty line
            a line starting with #, containing meta-info for a snapshot (strain, strain rate, etc)
            a bunch of lines describing a snapshot
            an empty line
            a line starting with #, containing meta-info for a snapshot (strain, strain rate, etc)
            a bunch of lines describing a snapshot
            ...
        Examples include par_ and int_ files.

    Parameters:
        in_file: the filename, or anything that can be taken as a first argument to pd.read_table
        field_nb: the number of fields (columns) in a snapshot. If not provided, the field nb is guessed from the file name.

    Returning values:
        frames: a list of snapshots
        strains_: the associated strains
        shear_rates_: the associated strain rates
    """
    field_nb_d = {'par': 15, 'int': 17}
    if field_nb is None:
        file_type = in_file[in_file.rfind('/')+1:in_file.find('_')]
        field_nb = field_nb_d[file_type]

    names = [str(i) for i in range(1,field_nb+1)]
    frames = pd.read_table(in_file, delim_whitespace=True, names=names, skiprows=field_nb+6)

    framebreaks = np.nonzero((frames['1'] == '#').as_matrix())[0] # locate empty lines
    frames = frames.as_matrix()
    frame_metadata = frames[framebreaks][:,1:].astype(np.float)
    shear_rates_ = frame_metadata[:,2]
    strains_ = frame_metadata[:,0]
    framebreaks = framebreaks[1:]
    frames = np.split(frames, framebreaks)

    for i in range(len(frames)):
        frames[i] = frames[i][1:]

    return frames, strains_, shear_rates_

def read_data_file(fname):
    """
    Purpose:
        Read any LF_DEM file that has a one-time-step-one-line structure, i.e:
            a header (x lines starting with #)
            an empty line
            a bunch of lines of data, one per time-step describing a snapshot

            Examples include data_ and st_ files.
    Parameters:
        fname: the filename, or anything that can be taken as a first argument to np.genfromtxt

    Returning values:
        A numpy array containing the data
    """
    dat = np.genfromtxt(fname)
    return dat
