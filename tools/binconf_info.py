#!/usr/bin/env python

import struct
import tty, sys, termios

def misuse():
    misuse_string = """
          Utilisation: 
           %s conf_file_in_binary_fmt
"""
    print(misuse_string % (sys.argv[0]))
    sys.exit(1)


if len(sys.argv) < 2 :
    misuse()

filename = sys.argv[1]

with open(filename, mode='rb') as f:
    conf = f.read()

uisize = 2
isize = 4
dsize = 8
loc=0
np = struct.unpack("i",conf[loc:isize])[0]
loc += isize
vf = struct.unpack("d",conf[loc:loc+dsize])[0]
loc += dsize
lx = struct.unpack("d",conf[loc:loc+dsize])[0]
loc += dsize
ly = struct.unpack("d",conf[loc:loc+dsize])[0]
loc += dsize
lz = struct.unpack("d",conf[loc:loc+dsize])[0]
loc += dsize
lees_x = struct.unpack("d",conf[loc:loc+dsize])[0]
loc += dsize
lees_y = struct.unpack("d",conf[loc:loc+dsize])[0]
loc += dsize

print("Particle number : ", np)
print("Volume/Area fraction : ", vf)
print("lx : ", lx)
print("ly : ", ly)
print("lz : ", lz)
print("shear displacement (x,y): ", lees_x, ",", lees_y, "\n")

print("Print full configuration (y/n)?")

fd = sys.stdin.fileno()
old_settings = termios.tcgetattr(fd)

try:
    tty.setraw(sys.stdin.fileno())
    ch = sys.stdin.read(1)
finally:
    termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)

if ch=="y":
    for i in range(np):
        x = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        y = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        z = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        r = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        print(x,y,z,r)


    nc = struct.unpack("=I",conf[loc:isize])[0]
    loc += isize

    for i in range(nc):
        p0 = struct.unpack("i",conf[loc:uisize])[0]
        loc += uisize
        p1 = truct.unpack("i",conf[loc:uisize])[0]
        loc += uisize
        dtx = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        dty = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        dtz = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        drx = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        dry = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        drz = struct.unpack("d",conf[loc:loc+dsize])[0]
        loc += dsize
        
        print(p0,p1,dtx,dty,dtz,drx,dry,drz)
