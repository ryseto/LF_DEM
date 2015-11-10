#!/usr/bin/env python
#
# yapgen.py, part of LF_DEM
#
# Purpose: yapgen.py is a script to generate a Yaplot visualization
#          of a LF_DEM simulation from the par_* and int_* files generated
#          by LF_DEM.
#
#  Romain Mari, 2015

import sys
import pandas as pd
import numpy as np

def read_lf_dem_snapshots(in_file, field_nb):
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

    Parameters:
        in_file: the filename, or anything that can be taken as a first argument to pd.read_table
        field_nb: the number of fields (columns) in a snapshot

    Returning values:
        frames: a list of snapshots
        strains_: the associated strains
        shear_rates_: the associated strain rates
    """
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
#    frames = [ frame[1:] for frame in frames ] # remove header
#    print frames
    return frames, strains_, shear_rates_

def read_data(posfile, intfile):
    """
    Purpose:
        Read a par_ and int_ file.

    Returning values:
        pos_frames: a list of snapshots from the par_ file
        int_frames: a list of snapshots from the int_ file
        strains: the associated strains
        shear_rates: the associated strain rates
    """
    field_nb = 15

    # load with pandas read_table
    pos_frames, strains, shear_rates = read_lf_dem_snapshots(posfile, field_nb)

    # load with pandas read_table
    field_nb = 17
    int_frames, strains, shear_rates = read_lf_dem_snapshots(intfile, field_nb)
#    print pos_frames[:3], int_frames[:3]
    return pos_frames, int_frames, strains, shear_rates

def hfill_array(cmd_array):
    """
        Purpose:
            From an array of strings with shape (N,m), with m<=7,
            get an array of strings with shape (N,7).
    """
    fill_nb = 7 - cmd_array.shape[1]
    height = cmd_array.shape[0]
    filler = np.tile([''], (height,fill_nb))
    return np.column_stack((cmd_array,filler))

def cmd(switch_type, switch_value):
    """
        Purpose:
            Get an array of strings for commands of type switch_type taking values
            switch_value.
    """
    sval = np.array(switch_value, dtype = np.str,ndmin=1)
    switch_cmd = np.empty(sval.shape[0], dtype = np.str)
    switch_cmd[:] = switch_type
    cmd = np.column_stack((switch_cmd, sval))
    return hfill_array(cmd)

def add_cmd(yap_array, switch_type, switch_value):
    """
        Purpose:
            Append to yap_array an array of strings for commands of type switch_type taking values
            switch_value.
    """
    return np.row_stack((yap_array, cmd(switch_type, switch_value)))

def layer_switch(value):
    return cmd('y', value)
def add_layer_switch(yap_array, value):
    return add_cmd(yap_array, 'y', value)

def color_switch(value):
    return cmd('@', value)
def add_color_switch(yap_array, value):
    return add_cmd(yap_array, '@', value)

def radius_switch(value):
    return cmd('r', value)
def add_radius_switch(yap_array, value):
    return add_cmd(yap_array, 'r', value)


def pair_cmd_and_switch(cmd, switch):
    """
    Purpose:
        Common use case: you want to change state (e.g. width) for every object.
        You can do that in an array-like fashion, generating cmd and switch arrays separately,
        and blending them afterwards. This is what this function is for.
    """
    return np.reshape(np.column_stack((switch,cmd)),(2*switch.shape[0], switch.shape[1]))



pos_fname = sys.argv[1]
forces_fname = pos_fname.replace("par_", "int_")
positions, forces, strains, shear_rates = read_data(pos_fname, forces_fname)

yap_filename = pos_fname.replace("par_", "y_")
yap_file = open(yap_filename,'wb')

nb_of_frames = len(strains)
i = 0
for f,p,strain,rate in zip(forces,positions, strains, shear_rates):
    # for each interaction: the particle indices
    part1 = f[:,0].astype(np.int)
    part2 = f[:,1].astype(np.int)

    # for each interaction: the particle positions
    r1 = p[part1,2:5].astype(np.float)
    r2 = p[part2,2:5].astype(np.float)

    # exclude interactions across the boundaries
    keep = np.linalg.norm(r2-r1,axis=1) < 4.
    r1 = r1[keep]
    r2 = r2[keep]
    f = f[keep]


    # display a line joining the center of interacting particles
    # with a thickness proportional to the normal force
    normal_forces = (f[:,7]+f[:,11]+f[:,15]).astype(np.float)
    r1r2 = np.hstack((r1,r2)).astype(np.str)
    yap_out = layer_switch(2)
    yap_out = add_color_switch(yap_out,4)

    force_factor = 0.00001 # to convert the force to a thickness. case-by-case.
    interaction_sticks = cmd('s', r1r2)
    interaction_widths = cmd('r', force_factor*np.abs(normal_forces)) # this
    interaction_out = pair_cmd_and_switch(interaction_sticks, interaction_widths)
    yap_out = np.row_stack((yap_out, interaction_out))

    # display a circle for every particle
    yap_out = add_layer_switch(yap_out, 3)
    yap_out = add_color_switch(yap_out,3)
    r = p[:,2:5].astype(np.float)
    rad = p[:,1].astype(np.float)
    particle_circle_positions = cmd('c', r)
    particle_circle_radius = cmd('r', rad)
    particle_out = pair_cmd_and_switch(particle_circle_positions, particle_circle_radius)
    yap_out = np.row_stack((yap_out, particle_out))

    yap_out = add_layer_switch(yap_out, 5)

    yap_out = np.row_stack((yap_out, ['t',str(10.),str(0.),str(10.),'strain='+str(strain),'','']))
    np.savetxt(yap_file, yap_out, fmt="%s "*7)
    yap_file.write("\n".encode('utf-8'))
    i += 1
    print("\rframe "+str(i)+"/"+str(nb_of_frames),end="",flush=True)
yap_file.close()
