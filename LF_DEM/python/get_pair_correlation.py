#!/opt/local/bin/python
#coding=utf-8

import sys, os
import math

self_path=os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.append(self_path+'/cython')
import LF_DEM_posfile_reading
import pair_correlation


def misuse():
    misuse_string = """
          Utilisation: 
           %s -r  [ r_bin_nb r_min r_max INPUT_FILE ]
           %s -c  [ r_bin_nb theta_bin_nb r_min r_max INPUT_FILE ]
           %s -s  [ r_bin_nb theta_bin_nb phi_bin_nb r_min r_max INPUT_FILE ]
"""
    print misuse_string % (sys.argv[0], sys.argv[0], sys.argv[0])
    exit(1)

def init():
    global r_bin_nb, theta_bin_nb, phi_bin_nb, r_min, r_max, mode
    args_nb = len(sys.argv)

    if args_nb < 2 :
        misuse()

    option = sys.argv[1]

    if option == "-s":
        if args_nb != 8:
            misuse()
        r_bin_nb=int(sys.argv[2])
        theta_bin_nb=int(sys.argv[3])
        phi_bin_nb=int(sys.argv[4])
        r_min=float(sys.argv[5])
        r_max=float(sys.argv[6])
        input_stream=open(str(sys.argv[7]),"r")
        mode="s"

    if option == "-c":
        if args_nb != 7:
            misuse()
        r_bin_nb=int(sys.argv[2])
        theta_bin_nb=int(sys.argv[3])
        r_min=float(sys.argv[4])
        r_max=float(sys.argv[5])
        input_stream=open(str(sys.argv[6]),"r")
        mode="c"


    if option == "-r":
        if args_nb != 6:
            misuse()
        r_bin_nb=int(sys.argv[2])
        r_min=float(sys.argv[3])
        r_max=float(sys.argv[4])
        input_stream=open(str(sys.argv[5]),"r")
        mode="r"

    return input_stream


mode=str()
r_bin_nb=int()
theta_bin_nb=int()
phi_bin_nb=int()
r_max=float()
r_min=float()
stream=init()

pos_stream=LF_DEM_posfile_reading.Pos_Stream(stream)

if pos_stream.dimension() == 2 and mode == "s":
    sys.stderr.write(' WARNING : input is 2d, switching to circular coordinates mode ')
    mode = "c"

if mode == "r":
    params = [r_bin_nb, r_min, r_max]
if mode == "c":
    params = [r_bin_nb,theta_bin_nb, r_min, r_max]
if mode == "s":
    params = [r_bin_nb,theta_bin_nb,phi_bin_nb, r_min, r_max]


twopoint_correl=pair_correlation.PairCorrelation( mode , params )
snapshot_nb=0

pos_stream.get_snapshot()
while pos_stream.get_snapshot():
    sys.stderr.write('Strain analyzed : '+str(pos_stream.time())+'\n')

    twopoint_correl.update_field(pos_stream)

twopoint_correl.normalize(pos_stream.np(), pos_stream.rho())
twopoint_correl.print_to(sys.stdout)


