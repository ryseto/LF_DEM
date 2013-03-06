cimport libc.math
import sys

def cart2sph(x):
    rtp=[]
    return cart2sph_nolist(x[0], x[1], x[2])


# r, theta \in [0,pi], phi \in [0, 2pi]
def cart2sph_nolist(double x, double y , double z): 
    rtp=[]

    cdef double theta
    cdef double phi
    cdef double coord
    cdef double pi=3.141592653589793238462643
    cdef double cosphi


    cdef double r=x*x+y*y+z*z
    r=libc.math.sqrt(r)
    rtp.append(r)
    theta=libc.math.acos(z/r)
    rtp.append(theta)

    if theta<1e-10: # deal with singularity
        phi=1e-10
    else:
        cosphi=x/r/libc.math.sin(theta)
        if libc.math.fabs(cosphi)<1.:
            phi=libc.math.acos(cosphi)
        else:
            if libc.math.fabs(cosphi)-1.<1e-10: # rounding error leading to |cos(phi)| = 1 + epsilon
                if x/libc.math.sin(theta)<0.:
                    phi=pi-1e-10
                else:
                    phi=1e-10
            else:
                sys.stderr.write(' cart2sph ERROR : could not compute phi. (x,y,z) is ('+str(x)+', '+str(y)+', '+str(z)+')'+' (r,theta) is ('+str(r)+', '+str(theta)+') \n')
                sys.exit(-1)
	    
    if y<0.:
        phi=2.*pi-phi  # phi \in [0, 2pi]


    rtp.append(phi)

    return rtp


def cart2circ(x):
    return cart2circ_nolist(x[0], x[1])


# r, theta \in [0,2*pi]
def cart2circ_nolist(double x, double y): 
    rt=[]

    cdef double theta
    cdef double coord
    cdef double pi=3.141592653589793238462643


    cdef double r=x*x+y*y
    r=libc.math.sqrt(r)
    rt.append(r)
    theta=libc.math.acos(x/r)

    if y<0.:
        theta=2.*pi-theta  # theta \in [0, 2pi]


    rt.append(theta)

    return rt


