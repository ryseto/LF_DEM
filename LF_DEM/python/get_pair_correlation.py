#!/usr/local/bin/python
#coding=utf-8

import sys
import math
sys.path.append('./cython')
import LF_DEM_posfile_reading
import pair_correlation



def init():
    global r_bin_nb, theta_bin_nb, phi_bin_nb, r_max, N, phi
    if len(sys.argv) != 8 :
        print "   Utilisation: ", sys.argv[0], "r_bin_nb theta_bin_nb phi_bin_nb r_max INPUT_FILE N phi"
        exit(1)

    r_bin_nb=int(sys.argv[1])
    theta_bin_nb=int(sys.argv[2])
    phi_bin_nb=int(sys.argv[3])
    r_max=float(sys.argv[4])
    input_stream=open(str(sys.argv[5]),"r")
    N=int(sys.argv[6])
    phi=float(sys.argv[7])


    return input_stream



r_bin_nb=int()
theta_bin_nb=int()
phi_bin_nb=int()
r_max=float()
N=int()
phi=float()
stream=init()

pos_stream=LF_DEM_posfile_reading.Pos_Stream(stream, N,phi)

twopoint_correl=pair_correlation.PairCorrelation(r_bin_nb,theta_bin_nb,phi_bin_nb, r_max)
snapshot_nb=0

pos_stream.get_snapshot()
while pos_stream.get_snapshot():
    sys.stderr.write(str(pos_stream.time())+'\n')

    twopoint_correl.update_field(pos_stream)

twopoint_correl.normalize()
twopoint_correl.print_to(sys.stdout)


