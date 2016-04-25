# lfdem_file.py, part of LF_DEM
#
#  Romain Mari, 2015
"""lfdem_file is a module providing functions
    to read output files from lf_dem"""

import pandas as pd
import numpy as np
import struct


def parse_file_metafield(key, value):
    if key == 'data':
        value = " ".join([key]+value)
        key = 'units'
    elif key == 'np':
        value = int(value[0])
    elif len(key) == 2 and key[0] == 'L':
        value = float(value[0])
    elif key == 'VF':
        value = float(value[0])
    elif key == 'LF_DEM':
        key = 'LF_DEM version'
        value = value[-1]
    return key, value


def parse_column_def(col_label_string, col_def):
    c = col_label_string
    c = c[c.find("#")+1:c.find(":")]
    col_def = " ".join(col_def)

    return col_def, c


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
    file_metadata['column def'] = {}

    header_len = 0

    while True:
        line = in_file.readline()
        if line[0] != '#':
            break

        data_list = line.split()

        if data_list[0] == '#':
            key, value = parse_file_metafield(data_list[1], data_list[2:])
            file_metadata[key] = value
        else:
            key, value = parse_column_def(data_list[0], data_list[1:])
            file_metadata['column def'][key] = value
        header_len += 1

    header_len += 1
    cols = convert_columndef_to_indices(file_metadata['column def'])
    col_nb = 0
    for k in cols:
        value = cols[k]
        try:
            if col_nb < value.stop:
                col_nb = value.stop
        except AttributeError:
            if col_nb < value + 1:
                col_nb = value + 1
    in_file.seek(0, 0)

    return file_metadata, col_nb, header_len


def convert_columndef_to_indices(columndef_dict):
    c = dict(columndef_dict)
    for k in c:
        c[k] = c[k].split("-")
        if len(c[k]) > 1:
            c[k] = slice(int(c[k][0])-1, int(c[k][1]))
        else:
            c[k] = int(c[k][0])-1
    return c


def __read_snapshot_file_no_framemeta(in_file, field_nb, usecols, header_len):
    # read col0 separately to determine frame breaks
    col0 = np.genfromtxt(in_file, comments=' ', skip_header=header_len)
    framebreaks_col0 = np.nonzero(np.isnan(col0))[0]

    col0 = col0[np.logical_not(np.isnan(col0))].astype(np.float)

    # now read cols>0
    in_file.seek(0, 0)
    # merge it with col0 if necessary
    if usecols != "all":
        if usecols[0] == 0:  # assume ordered
            cols = np.genfromtxt(in_file, skip_header=header_len,
                                 usecols=usecols[1:])
            cols = np.column_stack((col0, cols))
        else:
            cols = np.genfromtxt(in_file, skip_header=header_len,
                                 usecols=usecols)
    else:
        cols = np.genfromtxt(in_file, skip_header=header_len)
    # get the frame breaks
    framebreaks_cols = framebreaks_col0 - np.arange(len(framebreaks_col0))
    framebreaks_cols = framebreaks_cols[1:]  # 1st break is line 0
    # split the data as a list of frames
    cols = np.split(cols, framebreaks_cols)
    return cols


def __read_snapshot_file_with_framemeta(in_file, field_nb, header_len):
    names = [str(i) for i in range(1, field_nb+1)]
    frames = pd.read_table(in_file, delim_whitespace=True,
                           names=names, skiprows=header_len)

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


def read_snapshot_file(fname, usecols="all", frame_meta=True):
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
        usecols: which columns to read, optional (default all columns)
        frame_meta: get the metadata of each frame, optional (default True)
                    [Note that for large files getting the metadata
                     can be very memory consuming]
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

    file_metadata, field_nb, header_len = get_file_metadata(in_file)
    in_file.close()
    if frame_meta:
        in_file = open(fname, "r")
    else:
        in_file = open(fname, "rb")  # for genfromtxt

    if frame_meta:
        return __read_snapshot_file_with_framemeta(in_file,
                                                   field_nb,
                                                   header_len)\
                + (file_metadata,)
    else:
        return __read_snapshot_file_no_framemeta(in_file,
                                                 field_nb,
                                                 header_len,
                                                 usecols=usecols),\
                file_metadata


def read_data_file(fname, usecols="all"):
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
        usecols: which columns to read, optional (default all columns)

    Returning values:
        data: a numpy array containing the data
        metadata: the files metadata
    """
    metadata, col_nb = get_file_metadata(fname)
    if usecols == "all":
        data = np.genfromtxt(fname)
    else:
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

    dat = np.genfromtxt(fname, usecols=np.arange(4))
    pos = dat[:, :3].astype(np.float)
    rad = dat[:, 3].astype(np.float)
    if openedfile:
        in_file.close()
    return pos, rad, meta_data


def popValue(t, stream):
    size = struct.calcsize(t)
    buf = stream.read(size)
    return struct.unpack(t, buf)


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

    return config
