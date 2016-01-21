#!/usr/bin/env python
import numpy as np
import sys
import lfdem_file as lf

if len(sys.argv) != 2:
    print(sys.argv[0], " conf\n")

fname = sys.argv[1]
pos, rad, meta = lf.read_conf_file(fname)
# get the top and bottom layers
ztop = pos[:,2]+rad
top_particles = ztop > float(meta['lz'])

zbottom = pos[:,2]-rad
bottom_particles = zbottom < 0

# freeze them
fixed_part = np.logical_or(bottom_particles, top_particles)

# now trick LF_DEM PBCs along z by extending the system along z
# and putting particles in the middle of this extended box,
zextension = 4*np.amax(rad)
print(meta['lz'])
meta['lz'] = str(float(meta['lz']) + zextension)
pos[:,2] += zextension/2
print(meta['lz'])

# separate frozen from mobile particles
fixed_pos = pos[fixed_part]
fixed_rad = rad[fixed_part]
mobile_pos = pos[np.logical_not(fixed_part)]
mobile_rad = rad[np.logical_not(fixed_part)]
mobile = np.column_stack((mobile_pos, mobile_rad))
out_data = mobile
fixed = np.column_stack((fixed_pos, fixed_rad))
out_data = np.row_stack((out_data, fixed))

# output all that
out_fname = fname.replace(".dat", "walls.dat")
with open(out_fname, "w") as out:
    line1 = '# '+' '.join(meta.keys())+'\n'
    line2 = '# '+' '.join(meta.values())+'\n'
    out.write(line1)
    out.write(line2)

out = open(out_fname, "ab")
np.savetxt(out, out_data)
out.close()
