#!/usr/bin/env python
#
# yapgen.py, part of LF_DEM
#
# Purpose: yapgen.py is a script to generate a Yaplot visualization
#          of a LF_DEM simulation from the par_* and int_* files generated
#          by LF_DEM.
#
#  Romain Mari, 2015

from __future__ import print_function
import sys
import pandas as pd
import numpy as np
import lfdem_file as lf
import pyaplot as pyp

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
    pos_frames, strains, shear_rates = lf.read_snapshot_file(posfile, field_nb)

    # load with pandas read_table
    field_nb = 17
    int_frames, strains, shear_rates = lf.read_snapshot_file(intfile, field_nb)
#    print pos_frames[:3], int_frames[:3]
    return pos_frames, int_frames, strains, shear_rates


def snaps2yap(pos_fname, force_factor):
    forces_fname = pos_fname.replace("par_", "int_")
    positions, forces, strains, shear_rates = read_data(pos_fname, forces_fname)

    yap_filename = pos_fname.replace("par_", "y_")
    yap_file = open(yap_filename,'wb')

    nb_of_frames = len(strains)
    i = 0
    for f,p,strain,rate in zip(forces,positions, strains, shear_rates):
        r1r2 = pyp.get_interaction_end_points(f,p)
        f, r1r2 = pyp.filter_interactions_crossing_PBC(f,r1r2)

        # display a line joining the center of interacting particles
        # with a thickness proportional to the normal force
        normal_forces = (f[:,7]+f[:,11]+f[:,15]).astype(np.float)
        normal_forces = force_factor*np.abs(normal_forces) # to convert the force to a thickness. case-by-case.
        yap_out = pyp.get_interactions_yaparray(r1r2, normal_forces)

        # display a circle for every particle
        pos = p[:,2:5].astype(np.float)
        rad = p[:,1].astype(np.float)
        yap_out = np.row_stack((yap_out, pyp.get_particles_yaparray(pos, rad)))

        yap_out = pyp.add_layer_switch(yap_out, 5)

        # display strain
        yap_out = np.row_stack((yap_out, ['t',str(10.),str(0.),str(10.),'strain='+str(strain),'','']))

        # output
        np.savetxt(yap_file, yap_out, fmt="%s "*7)
        yap_file.write("\n".encode('utf-8'))
        i += 1
        try:
            print("\rframe "+str(i)+"/"+str(nb_of_frames),end="",flush=True)
        except TypeError: # to work with Python 2.7.* importing print_function without flush arg
            print("\rframe "+str(i)+"/"+str(nb_of_frames),end="")
    yap_file.close()

def conf2yap(conf_fname):
    yap_filename = pos_fname.replace(".dat", ".yap")
    yap_file = open(yap_filename,'wb')

    positions, radii, meta = lf.read_conf_file(conf_fname)
    positions[:,0] -= float(meta['lx'])/2
    positions[:,1] -= float(meta['ly'])/2
    positions[:,2] -= float(meta['lz'])/2

    yap_out = pyp.get_particles_yaparray(positions, radii)
    # print(yap_out)
    np.savetxt(yap_file, yap_out, fmt="%s "*7)
    yap_file.write("\n".encode('utf-8'))
    yap_file.close()

if len(sys.argv) < 2:
    print(sys.argv[0], " par_or_conf_file [force_factor]\n")
    exit(1)

pos_fname = sys.argv[1]

if pos_fname.find("par_") > -1:
    if len(sys.argv)>2:
        force_factor = float(sys.argv[2])
    else:
        force_factor = 0.01
    snaps2yap(pos_fname, force_factor)
else:
    conf2yap(pos_fname)
