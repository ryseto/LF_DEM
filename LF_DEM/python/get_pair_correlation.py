#!/opt/local/bin/python
#coding=utf-8

import sys
import math
sys.path.append('./cython')
import LF_DEM_posfile_reading
import pair_correlation



def init():
    global r_bin_nb, theta_bin_nb, phi_bin_nb, r_min, r_max
    if len(sys.argv) != 7 :
        print "   Utilisation: ", sys.argv[0], "r_bin_nb theta_bin_nb phi_bin_nb r_min r_max INPUT_FILE"
        exit(1)


    r_bin_nb=int(sys.argv[1])
    theta_bin_nb=int(sys.argv[2])
    phi_bin_nb=int(sys.argv[3])
    r_min=float(sys.argv[4])
    r_max=float(sys.argv[5])
    input_stream=open(str(sys.argv[6]),"r")

    return input_stream



r_bin_nb=int()
theta_bin_nb=int()
phi_bin_nb=int()
r_max=float()
r_min=float()
stream=init()

pos_stream=LF_DEM_posfile_reading.Pos_Stream(stream)

if pos_stream.dimension() == 2:
    params = [r_bin_nb,theta_bin_nb, r_min, r_max]
if pos_stream.dimension() == 3:
    params = [r_bin_nb,theta_bin_nb,phi_bin_nb, r_min, r_max]

twopoint_correl=pair_correlation.PairCorrelation( pos_stream.dimension(), params )
snapshot_nb=0

pos_stream.get_snapshot()
while pos_stream.get_snapshot():
    sys.stderr.write(str(pos_stream.time())+'\n')

    twopoint_correl.update_field(pos_stream)

twopoint_correl.normalize(pos_stream.np(), pos_stream.rho())
twopoint_correl.print_to(sys.stdout)


