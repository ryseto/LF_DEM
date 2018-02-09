#!/usr/bin/env python
import numpy as np
import sys
import clfdem_file as lf

if len(sys.argv) != 3:
    print(" Usage: ", sys.argv[0], " conf wall_thickness\n")
    exit(1)
fname = sys.argv[1]
wall_thickness = float(sys.argv[2])

pos, rad, meta = lf.read_conf_file(fname)
# get the top and bottom layers

ztop = pos[:, 2]+wall_thickness
top_particles = ztop > float(meta['lz'])

zbottom = pos[:, 2]-wall_thickness*rad
bottom_particles = zbottom < 0

# freeze them
fixed_part = np.logical_or(bottom_particles, top_particles)

# now trick LF_DEM PBCs along z by extending the system along z
# and putting particles in the middle of this extended box,
zextension = 4*np.amax(rad)
meta['lz'] = str(float(meta['lz']) + zextension)
pos[:, 2] += zextension/2

# separate frozen from mobile particles
fixed_rad = rad[fixed_part]
mobile_pos = pos[np.logical_not(fixed_part)]
mobile_rad = rad[np.logical_not(fixed_part)]
mobile = np.column_stack((mobile_pos, mobile_rad))
# out_data = mobile
fixed_pos = pos[fixed_part]
fixed_vel = np.zeros(pos.shape)
fixed_vel[bottom_particles] = np.array([-1, 0, 0])
fixed_vel[top_particles] = np.array([1, 0, 0])
fixed_vel = fixed_vel[fixed_part]
fixed = np.column_stack((fixed_pos, fixed_rad, fixed_vel))

# out_data = np.row_stack((out_data, fixed))

# output all that
out_fname = fname.replace(".dat", "walls.dat")
with open(out_fname, "w") as out:
    line1 = '# '+' '.join(meta.keys())+' np_fixed\n'
    line2 = '# '+' '.join(meta.values())+' '+str(len(fixed_pos))+'\n'
    # print(line2)
    out.write(line1)
    out.write(line2)

out = open(out_fname, "ab")
np.savetxt(out, mobile)
np.savetxt(out, fixed)

out.close()
