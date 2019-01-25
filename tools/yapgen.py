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
import clfdem_file as clff
import lfdem_utils as lfu
import dict_utils as du
import pyaplot as pyp


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
    normal_forces_loc =\
        [a[1] for a in du.matching(coldef_dict, ["normal", "force"])]
    total_force = np.zeros(len(interactions), dtype=np.float32)
    for f in normal_forces_loc:
        total_force += interactions[:, f].astype(np.float32)

    return total_force



def particles_yaparray(pos, rad, angles=None):
    """
        Get yaplot commands (as an aray of strings) to display circles
        for each particle defined in (pos, rad).
        pos and rad must contain positions and radii of particles.
    """

    if pos.shape[-1]==3:
        npos = np.asarray(pos)
    else:
        npos = np.column_stack((np.zeros(len(pos)), pos))[:, [1, 0, 2]]
    
    particle_circle_positions = pyp.cmd('c', npos)
    particle_circle_radius = pyp.cmd('r', rad)
    yap_out = pyp.pair_cmd_and_switch(particle_circle_positions,
                                      particle_circle_radius)
    if angles is None:
        return yap_out
    else:
        # add crosses in 2d
        u1 = np.zeros(npos.shape)
        u2 = np.zeros(npos.shape)
        u1[:, 0] = np.cos(angles)
        u1[:, 2] = -np.sin(angles)
        u1[:, [0, 2]] *= rad[:, np.newaxis]
        u2[:, 0] = -u1[:, 2]
        u2[:, 2] = u1[:, 0]
        depth_shift = np.array([0, -0.1, 0])
        crosses = pyp.cmd('l', np.row_stack((np.hstack((npos + u1 + depth_shift, npos - u1 + depth_shift)),
                                             np.hstack((npos + u2 + depth_shift, npos - u2 + depth_shift)))))
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
    if r1r2.shape[-1]==4:
        zeros = np.zeros(len(r1r2))
        r1r2 = np.column_stack((r1r2[:,0], zeros, r1r2[:,[1,2]], zeros, r1r2[:,3]))
    f, r1r2 = filter_interactions_crossing_PBC(int_snapshot, r1r2)

    # display a line joining the center of interacting particles
    # with a thickness proportional to the normal force
    normal_forces = get_normal_force(f, icols)
    #  convert the force to a thickness. case-by-case.
    if f_factor is None and len(normal_forces)>0:
        f_factor = 0.5/np.max(np.abs(normal_forces))
    normal_forces = f_factor*np.abs(normal_forces)
    avg_force = np.mean(np.abs(normal_forces))
    large_forces = normal_forces > f_chain_thresh * avg_force

    try:
        # contacts
        contact_state = f[:, du.matching_uniq(icols, 'contact state')[1]]
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
    except KeyError:
        keep = np.ones(len(f), dtype=np.bool)
        yap_out = pyp.layer_switch(layer_contacts)
        yap_out = pyp.add_color_switch(yap_out, color_contacts)
        bonds = pyp.sticks_yaparray(r1r2[keep], normal_forces[keep])
        yap_out = np.row_stack((yap_out, bonds))

    return yap_out, f_factor


def cuboid(lx2, ly2, lz2):
    return np.array([[lx2, ly2, lz2, lx2, ly2, -lz2],
                     [lx2, ly2, lz2, lx2, -ly2, lz2],
                     [lx2, ly2, lz2, -lx2, ly2, lz2],
                     [lx2, -ly2, -lz2, lx2, -ly2, lz2],
                     [lx2, -ly2, -lz2, lx2, ly2, -lz2],
                     [lx2, -ly2, -lz2, -lx2, -ly2, -lz2],
                     [-lx2, ly2, -lz2, -lx2, ly2, lz2],
                     [-lx2, ly2, -lz2, -lx2, -ly2, -lz2],
                     [-lx2, ly2, -lz2, lx2, ly2, -lz2],
                     [-lx2, -ly2, lz2, -lx2, -ly2, -lz2],
                     [-lx2, -ly2, lz2, -lx2, ly2, lz2],
                     [-lx2, -ly2, lz2, lx2, -ly2, lz2]])


def rectangle(lx2, lz2):
    return np.array([[lx2, 0, lz2, lx2, 0, -lz2],
                     [lx2, 0, lz2, -lx2, 0, lz2],
                     [-lx2, 0, -lz2, -lx2, 0, lz2],
                     [-lx2, 0, -lz2, lx2, 0, -lz2]])


def snaps2yap(pos_fname,
              yap_file,
              f_factor=None,
              f_chain_thresh=None):

    forces_fname = pos_fname.replace("par_", "int_")
    par_f = clff.snapshot_file(pos_fname)
    int_f = clff.snapshot_file(forces_fname)


    pcols = par_f.column_def()
    icols = int_f.column_def()

    is2d = float(par_f.meta_data()['Ly']) == 0
    i = 0
    frame_int = int_f.__next__()
    strain_int = du.matching_uniq(frame_int[0], ["cu".encode('utf8'), "strain".encode('utf8')])
    for frame_par in par_f:
        # *cu*rvilinear or *cu*mulated strain depending on LF_DEM version
        strain = du.matching_uniq(frame_par[0], ["cu".encode('utf8'), "strain".encode('utf8')])

        # display a circle for every particle
        yap_out = pyp.layer_switch(3)
        yap_out = pyp.add_color_switch(yap_out, 3)

        if 'position (x, y, z)' in pcols:
            pos_slice = pcols['position (x, y, z)']
        else:
            pos_slice = slice(pcols['position x'], pcols['position z']+1)
        pos = frame_par[1][:, pos_slice].astype(np.float)
        rad = frame_par[1][:, pcols['radius']].astype(np.float)
        if is2d:
            try:
                angle = frame_par[1][:, pcols['angle']].astype(np.float)
                particles, crosses = particles_yaparray(pos, rad, angles=angle)
                yap_out = np.row_stack((yap_out, particles))
            except KeyError:
                yap_out = np.row_stack((yap_out, particles_yaparray(pos, rad)))
        else:
            yap_out = np.row_stack((yap_out, particles_yaparray(pos, rad)))

        # display bounding box
        yap_out = pyp.add_layer_switch(yap_out, 4)
        yap_out = pyp.add_color_switch(yap_out, 0)
        lx2 = float(par_f.meta_data()['Lx'])/2
        ly2 = float(par_f.meta_data()['Ly'])/2
        lz2 = float(par_f.meta_data()['Lz'])/2
        if not is2d:
            yap_out = pyp.add_cmd(yap_out, 'l', cuboid(lx2, ly2, lz2))
        else:
            yap_out = pyp.add_cmd(yap_out, 'l', rectangle(lx2, lz2))

        # display strain
        yap_out = pyp.add_layer_switch(yap_out, 5)
        yap_out = pyp.add_color_switch(yap_out, 1)
        yap_out = np.row_stack((yap_out,
                                ['t', 1.1*lx2, str(0.), 1.1*lz2,
                                    'strain='+str(strain[1]), '', '']))

        # display interactions, if any
        while strain_int < strain:
            try:
                frame_int = int_f.__next__()
                strain_int = du.matching_uniq(frame_int[0], ["cu".encode('utf8'), "strain".encode('utf8')])
            except StopIteration:
                break
        if strain_int == strain:
            interaction_bonds, f_factor =\
            interactions_bonds_yaparray(frame_int[1], frame_par[1],
                                        icols, pcols,
                                        f_factor=f_factor,
                                        f_chain_thresh=f_chain_thresh,
                                        layer_contacts=1,
                                        layer_noncontacts=2,
                                        color_contacts=4,
                                        color_noncontacts=5)
            yap_out = np.row_stack((yap_out, interaction_bonds))


        # output
        np.savetxt(yap_file, yap_out, fmt="%s "*7)
        yap_file.write("\n".encode('utf-8'))
        i += 1
        out_str = "\r frame " + str(i) +\
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
    positions, radii, meta = clff.read_conf_file(conf_fname)
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
