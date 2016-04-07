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

    header_len = 0
    while True:
        line = in_file.readline()
        if line[0] != '#':
            break

        data_list = line.split()

        if data_list[0] == '#':
            file_metadata[data_list[1]] = data_list[2:]
        else:
            file_metadata['column def'] += line
        header_len += 1

    in_file.seek(0, 0)

    return file_metadata, header_len


def __read_snapshot_file_no_framemeta(in_file, field_nb, usecols):
    # read col0 separately to determine frame breaks
    header_len = field_nb+6
    col0 = np.genfromtxt(in_file, comments=' ', skip_header=header_len)
    framebreaks_col0 = np.nonzero(np.isnan(col0))[0]

    col0 = col0[np.logical_not(np.isnan(col0))].astype(np.float)

    # now read cols>0
    in_file.seek(0, 0)
    # merge it with col0 if necessary
    if usecols is not None and usecols[0] == 0:  # assume ordered
        cols = np.genfromtxt(in_file, skip_header=header_len,
                             usecols=usecols[1:])
        cols = np.column_stack((col0, cols))
    else:
        cols = np.genfromtxt(in_file, skip_header=header_len,
                             usecols=usecols)
    # get the frame breaks
    framebreaks_cols = framebreaks_col0 - np.arange(len(framebreaks_col0))
    framebreaks_cols = framebreaks_cols[1:]  # 1st break is line 0
    # split the data as a list of frames
    cols = np.split(cols, framebreaks_cols)
    return cols


def __read_snapshot_file_with_framemeta(in_file, field_nb, usecols):
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

    return (frames, strains_, shear_rates_, frame_metadata)


def read_snapshot_file(fname, field_nb=None, usecols=None, frame_meta=True):
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
        if frame_meta == False:
            frames: a list of snapshots
            file_metadata: the metadata of the whole file
        else:
            frames: a list of snapshots
            strains_: the associated strains
            shear_rates_: the associated strain rates
            frame_metadata: the metadata for each snapshot
            file_metadata: the metadata of the whole file
    """

    try:
        in_file = open(fname, "r")
    except TypeError:
        in_file = fname

    file_metadata, header_len = get_file_metadata(in_file)
    in_file.close()
    if frame_meta:
        in_file = open(fname, "r")
    else:
        in_file = open(fname, "rb")  # for genfromtxt

    if field_nb is None:
        field_nb = 0
        field_nb += header_len - 6

    if frame_meta:
        return __read_snapshot_file_with_framemeta(in_file,
                                                   field_nb,
                                                   usecols)\
                + (file_metadata,)
    else:
        return __read_snapshot_file_no_framemeta(in_file,
                                                 field_nb,
                                                 usecols=usecols),\
                file_metadata


def read_data_file(fname, usecols=None):
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
    metadata, header_len = get_file_metadata(fname)
    data = np.genfromtxt(fname, usecols=usecols)
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
