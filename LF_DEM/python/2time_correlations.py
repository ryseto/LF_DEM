#!/opt/local/bin/python
#coding=utf-8

# compile with
# cython *.pyx
# gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.7 -Wl,-Bsymbolic-functions -Wl,-z,relro
#

import sys
import math
sys.path.append('./cython')
import LF_DEM_posfile_reading
import two_time_correlator as Cor


def init():
    if not visual:
        arg_nb=3
        if len(sys.argv) != arg_nb :
            print "   Utilisation: ", sys.argv[0], "INPUT_FILE OUTPUT_files"
            exit(1)
            
        input_stream=open(str(sys.argv[1]),"r")
        out_filenames=sys.argv[2]
        
        return [ input_stream, out_filenames]
    else:
        arg_nb=2
        if len(sys.argv) != arg_nb :
            print "   Utilisation: ", sys.argv[0], "INPUT_FILE"
            exit(1)
            
        input_stream=open(str(sys.argv[1]),"r")
        
        return input_stream
        

def print_out():
    
    imported_cor = []
    for cor in correlators:
        imported_cor.append(cor.export_correlator())


    out_stream=open(out_fname, "w")
    title = """
    Two-time correlations + standard deviation ("X4"). 
    column #1 is "time" difference, in snapshot numbers.
    There are %i wave_vectors stored here.
    These are:
    """
    col_descriptor = " column #%i and #%i, k_x = %f, k_y = %f, k_z = %f \n"

    out_stream.write(title)
    col_nb=2
    for cor in correlators:
        wv = cor.wave_vector()
        out_stream.write(col_descriptor % (col_nb, col_nb+1, wv[0], wv[1], wv[2])  )
        col_nb += 2

    for point in sorted(imported_cor[0].iterkeys()):
        out_stream.write(str(point)+' ')
        for cor in imported_cor:
            out_stream.write(str(cor[point][0])+' '+str(cor[point][1])+' ')

        out_stream.write('\n')
    out_stream.close()

    


#================= Main ========================#

# two modes:
# - if visual is False, computes 2 time correlators and stores it
#   in files which name are given in input
# - if visual is True, just outputs positions minus advection
visual=False

# print time currently treated
verbose=False

# Parameters for visual=True :
# update_period : point sampling for correlators, in snapshot unit
# starting_period : for averaging, CorrelatorSet starts a new correlator every starting_period, in snapshot unit
# max_time : maximum time difference to be considered in correlators, in snapshot unit
# output_period : ouput averaged correlators every output_period, in snapshot unit
update_period = 2
starting_period = 100
max_time = 100
output_period = 100

if not visual:
    [ stream, out_fname ]=init()
else:
    stream=init()



pos_stream=LF_DEM_posfile_reading.Pos_Stream(stream)


# generate wave vectors in shear plane
# imax -> nb of wave vectors
# amp -> their norm
# store them in wave_vectors

wave_vectors=[]
imax=6
pi=3.141592653589793238462643
correlators=[]
amp=6.5/1.
for i in range(imax):
    theta=float(i)*pi/float(imax)
    wave_vectors.append([amp*math.cos(theta),0.,amp*math.sin(theta)])




pos_stream.get_snapshot()
positions=pos_stream.positions_deepcopy()

previous_update_time=0
previous_output_time=0

for k in wave_vectors: 
    correlators.append(Cor.TwoTimeCorrelatorSet(k, starting_period, max_time))

if visual:
    for i in positions:
        for s in range(3):
            positions[i][s]=0.

while pos_stream.get_snapshot():

    for i in positions:
        disp=pos_stream.non_affine_displacement(i)
        for s in range(3):
            positions[i][s]+=disp[s]

    if visual:
        sys.stderr.write(str(pos_stream.time())+'\n')
        print ""
        print "                         realtime = ", pos_stream.time()
        for i in positions:
            print i, positions[i][0], positions[i][1], positions[i][2],' 0 0 0'

    if not visual:
        if pos_stream.time_index() - previous_update_time == update_period:
            for cor in correlators:
                cor.update(positions, pos_stream.time_index())
            previous_update_time=pos_stream.time_index()    
            if verbose:
                sys.stderr.write(' time : '+str(pos_stream.time())+', time index : '+str(pos_stream.time_index())+'\n')

        if pos_stream.time_index() - previous_output_time == output_period:
            previous_output_time = pos_stream.time_index()
            print_out()
