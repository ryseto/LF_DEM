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
import lfdem_utils as lfu
import dict_utils as du
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


def filter_interactions_crossing_PBC(f, r1r2, cutoff=4):
    """
        Exclude interactions across the boundaries.
        Return values of f and r1r2 where norm(r1r2[:,3:]-r1r2[:,:3])<cutoff.
    """
    r1 = r1r2[..., :3]
    r2 = r1r2[..., 3:]
    keep = np.linalg.norm(r2-r1, axis=-1) < cutoff
    r1r2 = r1r2[keep]
    f = f[keep]
    return f, r1r2


def get_normal_force(interactions, coldef_dict):
    lub_loc = coldef_dict['normal part of the lubrication force']
    lub_force = interactions[:, lub_loc]
    cont_loc = coldef_dict['norm of the normal part of the contact force']
    contact_force = interactions[:, cont_loc]
    rep_loc = coldef_dict['norm of the normal repulsive force']
    repulsive_force = interactions[:, rep_loc]
    return (lub_force + contact_force + repulsive_force).astype(np.float32)


def particles_yaparray(pos, rad, angles=None):
    """
        Get yaplot commands (as an aray of strings) to display circles
        for each particle defined in (pos, rad).
        pos and rad must contain positions and radii of particles.
    """

    particle_circle_positions = pyp.cmd('c', pos)
    particle_circle_radius = pyp.cmd('r', rad)
    yap_out = pyp.pair_cmd_and_switch(particle_circle_positions,
                                      particle_circle_radius)
    if angles is None:
        return yap_out
    else:
        # add crosses in 2d
        u1 = -np.ones(pos.shape)   # so that they appear in front
        u2 = -np.ones(pos.shape)
        u1[:, 0] = np.cos(angles)
        u1[:, 2] = np.sin(angles)
        u1[:, [0, 2]] *= rad[:, np.newaxis]
        u2[:, 0] = -u1[:, 2]
        u2[:, 2] = u1[:, 0]
        crosses = pyp.cmd('l', np.row_stack((np.hstack((pos+u1, pos-u1)),
                                             np.hstack((pos+u2, pos-u2)))))
        return yap_out, crosses


def interactions_bonds_yaparray(int_snapshot,
                                par_snapshot,
                                icols,
                                pcols,
                                f_factor=None,
                                f_chain_thresh=None,
                                layer_contacts=1,
                                layer_noncontacts=2,
                                color_contacts=1,
                                color_noncontacts=2):
    if f_chain_thresh is None:
        f_chain_thresh = 0

    r1r2 = lfu.get_interaction_end_points(int_snapshot, par_snapshot,
                                          icols, pcols)
    f, r1r2 = filter_interactions_crossing_PBC(int_snapshot, r1r2)

    # display a line joining the center of interacting particles
    # with a thickness proportional to the normal force
    normal_forces = get_normal_force(f, icols)
    #  convert the force to a thickness. case-by-case.
    if f_factor is None:
        f_factor = 0.5/np.max(np.abs(normal_forces))
    normal_forces = f_factor*np.abs(normal_forces)
    avg_force = np.mean(np.abs(normal_forces))
    large_forces = normal_forces > f_chain_thresh * avg_force

    contact_state = f[:, du.matching_uniq(icols, 'contact state')[1]]

    # contacts
    keep = np.logical_and(contact_state > 0, large_forces)
    yap_out = pyp.layer_switch(layer_contacts)
    yap_out = pyp.add_color_switch(yap_out, color_contacts)
    contact_bonds = pyp.sticks_yaparray(r1r2[keep], normal_forces[keep])
    yap_out = np.row_stack((yap_out, contact_bonds))

    # non contacts
    keep = np.logical_and(contact_state == 0, large_forces)
    yap_out = pyp.add_layer_switch(yap_out, layer_noncontacts)
    yap_out = pyp.add_color_switch(yap_out, color_noncontacts)
    non_contact_bonds = pyp.sticks_yaparray(r1r2[keep], normal_forces[keep])
    yap_out = np.row_stack((yap_out, non_contact_bonds))

    return yap_out, f_factor


def snaps2yap(pos_fname,
              yap_file,
              f_factor=None,
              f_chain_thresh=None):

    forces_fname = pos_fname.replace("par_", "int_")
    positions, forces, strains, shear_rates, meta_pos, meta_int =\
        read_data(pos_fname, forces_fname)
    pcols = lf.convert_columndef_to_indices(meta_pos['column def'])
    icols = lf.convert_columndef_to_indices(meta_int['column def'])

    nb_of_frames = len(strains)
    is2d = meta_pos['Ly'] == 0
    i = 0
    for f, p, strain, rate in zip(forces, positions, strains, shear_rates):
        yap_out, f_factor =\
            interactions_bonds_yaparray(f, p, icols, pcols,
                                        f_factor=f_factor,
                                        f_chain_thresh=f_chain_thresh,
                                        layer_contacts=1,
                                        layer_noncontacts=2,
                                        color_contacts=4,
                                        color_noncontacts=5)

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
            particles, crosses = particles_yaparray(pos, rad, angles=angle)
            yap_out = np.row_stack((yap_out, particles))
            yap_out = pyp.add_color_switch(yap_out, 1)
            yap_out = np.row_stack((yap_out, crosses))
        else:
            yap_out = np.row_stack((yap_out, particles_yaparray(pos, rad)))

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
                  "   [force factor " + str(f_factor) + "]"
        try:
            print(out_str, end="", flush=True)
        except TypeError:
            # to work with Python 2.7.* importing print_function without flush
            print(out_str, end="")
    print("")


def conf2yap(conf_fname, yap_filename):
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
                                particles_yaparray(pos_mobile, rad_mobile)))
        yap_out = pyp.add_layer_switch(yap_out, 4)
        yap_out = pyp.add_color_switch(yap_out, 4)
        yap_out = np.row_stack((yap_out,
                                particles_yaparray(pos_fixed, rad_fixed)))
    else:
        yap_out = pyp.layer_switch(3)
        yap_out = pyp.add_color_switch(yap_out, 3)
        yap_out = np.row_stack((yap_out,
                                particles_yaparray(positions, radii)))

    pyp.savetxt(yap_filename, yap_out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-ff', '--force-factor', type=float)
    parser.add_argument('-ft', '--force-threshold', type=float)
    parser.add_argument('-o', '--output')
    parser.add_argument('file')

    args = vars(parser.parse_args(sys.argv[1:]))

    if args['file'].find("par_") > -1:
        if args['output'] is None:
            yap_filename =\
                args['file'].replace("par_", "y_").replace(".dat", ".yap")
        else:
            yap_filename = args['output']
        yap_file = open(yap_filename, 'wb')
        snaps2yap(args['file'],
                  yap_file,
                  f_factor=args['force_factor'],
                  f_chain_thresh=args['force_threshold'])
        yap_file.close()

    else:
        if args['output'] is None:
            yap_filename = args['file'].replace(".dat", ".yap")
        else:
            yap_filename = args['output']
        conf2yap(args['file'], yap_filename)
