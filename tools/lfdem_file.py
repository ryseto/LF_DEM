# lfdem_file.py, part of LF_DEM
#
#  Romain Mari, 2015
"""lfdem_file is a module providing functions
    to read output files from lf_dem"""

import pandas as pd
import numpy as np


def get_file_metadata(fname):
    """
    Purpose:
        Get the metadata of a file (the data contained in the header,
        lines starting with '#').
    """
    try:
        in_file = open(fname, "r")
    except TypeError:
        in_file = fname

    file_metadata = {}
    file_metadata['column def'] = str()

    while True:
        line = in_file.readline()
        if line[0] != '#':
            break

        data_list = line.split()

        if data_list[0] == '#':
            file_metadata[data_list[1]] = data_list[2:]
        else:
            file_metadata['column def'] += line

    in_file.seek(0, 0)

    return file_metadata


def read_snapshot_file(fname, field_nb=None):
    """
    Purpose:
        Read any LF_DEM file that has a "snapshot" structure, i.e. made of
        the following:
            * a header (x lines starting with #)
            * an empty line
            * a line starting with #, containing meta-info for a snapshot
                                                  (strain, strain rate, etc)
            * a bunch of lines describing a snapshot
            * an empty line
            * a line starting with #, containing meta-info for a snapshot
                                                  (strain, strain rate, etc)
            * a bunch of lines describing a snapshot
            ...
        Examples include par_ and int_ files.

    Parameters:
        fname: the filename, or a file like object
        field_nb: the number of fields (columns) in a snapshot.
                  If not provided, the field nb is guessed from the file name.

    Returning values:
        frames: a list of snapshots
        strains_: the associated strains
        shear_rates_: the associated strain rates
        frame_metadata: the metadata for each snapshot
        file_metadata: the metadata of the whole file
    """

    openedfile = True
    try:
        in_file = open(fname, "r")
    except TypeError:
        in_file = fname
        openedfile = False

    file_metadata = get_file_metadata(in_file)

    if field_nb is None:
        if openedfile:
            field_nb_d = {'par': 15, 'int': 17}
            file_type = fname[fname.rfind('/')+1:fname.find('_')]
            field_nb = field_nb_d[file_type]
        else:
            field_nb = 15
            print("Warning, using default field_nb=15")
            
    names = [str(i) for i in range(1, field_nb+1)]
    frames = pd.read_table(in_file, delim_whitespace=True,
                           names=names, skiprows=field_nb+6)

    # locate empty lines
    framebreaks = np.nonzero((frames['1'] == '#').as_matrix())[0]
    frames = frames.as_matrix()

    frame_metadata = frames[framebreaks][:, 1:].astype(np.float)

    shear_rates_ = frame_metadata[:, 2]
    strains_ = frame_metadata[:, 0]
    framebreaks = framebreaks[1:]
    frames = np.split(frames, framebreaks)

    for i in range(len(frames)):
        frames[i] = frames[i][1:]

    return frames, strains_, shear_rates_, frame_metadata, file_metadata


def read_data_file(fname):
    """
    Purpose:
        Read any LF_DEM file that has a one-time-step-one-line structure, i.e:
            a header (x lines starting with #)
            an empty line
            a bunch of lines of data, one per time-step describing a snapshot

            Examples include data_ and st_ files.
    Parameters:
        fname: the filename, or anything that can be taken as a first argument
               to np.genfromtxt

    Returning values:
        data: a numpy array containing the data
        metadata: the files metadata
    """
    metadata = get_file_metadata(fname)
    data = np.genfromtxt(fname)
    return data, metadata


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

    dat = np.genfromtxt(fname)
    pos = dat[:, :3].astype(np.float)
    rad = dat[:, 3].astype(np.float)
    if openedfile:
        in_file.close()
    return pos, rad, meta_data
