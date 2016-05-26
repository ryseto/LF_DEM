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
import numpy as np
import argparse
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
    # load with pandas read_table
    pos_frames, strains, shear_rates, dumb, meta_pos =\
        lf.read_snapshot_file(posfile)

    # load with pandas read_table
    int_frames, strains, shear_rates, dumb, meta_int =\
        lf.read_snapshot_file(intfile)
#    print pos_frames[:3], int_frames[:3]
    return pos_frames, int_frames, strains, shear_rates, meta_pos, meta_int


def get_normal_force(interactions, coldef_dict):
    lub_loc = coldef_dict['normal part of the lubrication force']
    lub_force = interactions[:, lub_loc]
    cont_loc = coldef_dict['norm of the normal part of the contact force']
    contact_force = interactions[:, cont_loc]
    rep_loc = coldef_dict['norm of the normal repulsive force']
    repulsive_force = interactions[:, rep_loc]
    return (lub_force + contact_force + repulsive_force).astype(np.float32)


def snaps2yap(pos_fname, force_factor=None, force_chain_threshold=None):
    if force_chain_threshold is None:
        force_chain_threshold = 0

    forces_fname = pos_fname.replace("par_", "int_")
    positions, forces, strains, shear_rates, meta_pos, meta_int =\
        read_data(pos_fname, forces_fname)
    pcols = lf.convert_columndef_to_indices(meta_pos['column def'])
    icols = lf.convert_columndef_to_indices(meta_int['column def'])

    yap_filename = pos_fname.replace("par_", "y_").replace(".dat", ".yap")
    yap_file = open(yap_filename, 'wb')

    nb_of_frames = len(strains)
    is2d = meta_pos['Ly'] == 0
    i = 0
    for f, p, strain, rate in zip(forces, positions, strains, shear_rates):
        r1r2 = pyp.get_interaction_end_points(f, p)
        f, r1r2 = pyp.filter_interactions_crossing_PBC(f, r1r2)

        # display a line joining the center of interacting particles
        # with a thickness proportional to the normal force
        normal_forces = get_normal_force(f, icols)
        #  convert the force to a thickness. case-by-case.
        if force_factor is None:
            force_factor = 0.5/np.max(np.abs(normal_forces))
        normal_forces = force_factor*np.abs(normal_forces)
        avg_force = np.mean(np.abs(normal_forces))
        large_forces = normal_forces > force_chain_threshold * avg_force

        contact_state = f[:, lf.strdict_get(icols, 'contact state')[1]]

        # contacts
        keep = np.logical_and(contact_state > 0, large_forces)
        yap_out = pyp.layer_switch(1)
        yap_out = pyp.add_color_switch(yap_out, 4)
        contact_bonds =\
            pyp.get_interactions_yaparray(r1r2[keep],
                                          normal_forces[keep])
        yap_out = np.row_stack((yap_out, contact_bonds))

        # non contacts
        keep = np.logical_and(contact_state == 0, large_forces)
        yap_out = pyp.add_layer_switch(yap_out, 2)
        yap_out = pyp.add_color_switch(yap_out, 5)
        non_contact_bonds =\
            pyp.get_interactions_yaparray(r1r2[keep],
                                          normal_forces[keep])
        yap_out = np.row_stack((yap_out, non_contact_bonds))

        # display a circle for every particle
        yap_out = pyp.add_layer_switch(yap_out, 3)
        yap_out = pyp.add_color_switch(yap_out, 3)
        if 'position (x, y, z)' in pcols:
            pos = p[:, pcols['position (x, y, z)']].astype(np.float)
        else:
            pos = p[:, pcols['position x']:pcols['position z']+1]\
                    .astype(np.float)
        rad = p[:, pcols['radius']].astype(np.float)
        if is2d:
            angle = p[:, pcols['angle']].astype(np.float)
            particles, crosses = \
                pyp.get_particles_yaparray(pos, rad, angles=angle)
            yap_out = np.row_stack((yap_out, particles))
            yap_out = pyp.add_color_switch(yap_out, 1)
            yap_out = np.row_stack((yap_out, crosses))
        else:
            yap_out = np.row_stack((yap_out,
                                    pyp.get_particles_yaparray(pos, rad)))

        yap_out = pyp.add_layer_switch(yap_out, 5)

        # display strain
        yap_out = np.row_stack((yap_out,
                                ['t', str(10.), str(0.), str(10.),
                                    'strain='+str(strain), '', '']))

        # output
        np.savetxt(yap_file, yap_out, fmt="%s "*7)
        yap_file.write("\n".encode('utf-8'))
        i += 1
        out_str = "\r frame " + str(i) + "/"+str(nb_of_frames) +\
                  " - " + yap_filename + \
                  "   [force factor " + str(force_factor) + "]"
        try:
            print(out_str, end="", flush=True)
        except TypeError:
            # to work with Python 2.7.* importing print_function without flush
            print(out_str, end="")
    print("")
    yap_file.close()


def conf2yap(conf_fname):
    yap_filename = pos_fname.replace(".dat", ".yap")
    print("Yap file : ", yap_filename)
    positions, radii, meta = lf.read_conf_file(conf_fname)
    positions[:, 0] -= float(meta['lx'])/2
    positions[:, 1] -= float(meta['ly'])/2
    positions[:, 2] -= float(meta['lz'])/2

    if 'np_fixed' in meta:
        # for conf with fixed particles
        split_line = len(positions) - int(meta['np_fixed'])
        pos_mobile, pos_fixed = np.split(positions, [split_line])
        rad_mobile, rad_fixed = np.split(radii, [split_line])
        yap_out = pyp.layer_switch(3)
        yap_out = pyp.add_color_switch(yap_out, 3)
        yap_out = np.row_stack((yap_out,
                                pyp.get_particles_yaparray(pos_mobile,
                                                           rad_mobile)))
        yap_out = pyp.add_layer_switch(yap_out, 4)
        yap_out = pyp.add_color_switch(yap_out, 4)
        yap_out = np.row_stack((yap_out,
                                pyp.get_particles_yaparray(pos_fixed,
                                                           rad_fixed)))
    else:
        yap_out = pyp.layer_switch(3)
        yap_out = pyp.add_color_switch(yap_out, 3)
        yap_out = np.row_stack(
                    (yap_out, pyp.get_particles_yaparray(positions, radii)))

    pyp.savetxt(yap_filename, yap_out)

parser = argparse.ArgumentParser()
parser.add_argument('-ff', '--force-factor', type=float)
parser.add_argument('-ft', '--force-threshold', type=float)
parser.add_argument('file')

args = vars(parser.parse_args(sys.argv[1:]))

if args['file'].find("par_") > -1:
    snaps2yap(args['file'],
              force_factor=args['force_factor'],
              force_chain_threshold=args['force_threshold'])
else:
    conf2yap(args['file'])
