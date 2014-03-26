cimport libc.math
from string import *
import sys
import copy
import numpy as np

cdef class Pos_Stream:

    cdef:
        instream
        is_file
        double lx, ly, lz, N, phi, V
        positions, old_positions, radius
        time_labels, last_read_time_label, current_time_index, first_time
        int dim
        pair_sep, pair_dist
		
    def __init__(self, input_stream):

        self.instream=input_stream
        try:
            self.instream.tell()
            self.is_file=True
        except IOError:
            self.is_file=False # means that we deal with stdin

        
        err_str = ' ERROR LF_DEM_posfile_reading : incorrect input file, '
        input_line=self.instream.readline()
        # get N
        input_line=self.instream.readline()
        fields=split(input_line)
        if fields[1] != 'np':
            sys.stderr.write(err_str+'no particle number \n')
        else:
            self.N = int(fields[2])

        # get VF
        input_line=self.instream.readline()
        fields=split(input_line)
        if fields[1] != 'VF':
            sys.stderr.write(err_str+'no volume fraction \n')
        else:
            self.phi = float(fields[2])

        # get Lx
        input_line=self.instream.readline()
        fields=split(input_line)
        if fields[1] != 'Lx':
            sys.stderr.write(err_str+'no Lx \n')
        else:
            self.lx = float(fields[2])

        # get Ly
        input_line=self.instream.readline()
        fields=split(input_line)
        if fields[1] != 'Ly':
            sys.stderr.write(err_str+'no Ly \n')
        else:
            self.ly = float(fields[2])

        # get Lz
        input_line=self.instream.readline()
        fields=split(input_line)
        if fields[1] != 'Lz':
            sys.stderr.write(err_str+'no Lz \n')
        else:
            self.lz = float(fields[2])

        if self.Ly() == 0.: # 2d case
            self.V = self.Lx()*self.Lz()
            self.dim=2
        else:               # 3d case
            self.V = self.Lx()*self.Ly()*self.Lz()
            self.dim=3

        self.reset_members()

    def __iter__(self):
        return self.positions.__iter__()

    def reset_members(self):
        self.positions=dict()
        self.old_positions=dict()
        self.radius=dict()
        self.last_read_time_label=[]
        self.time_labels=[]
        self.current_time_index=0
        test_input=self.get_snapshot()
        self.first_time=self.time()
        if test_input==False:
            print "Invalid input "
            sys.exit(1)

    def positions_copy(self):
        return self.positions.copy()

    def positions_deepcopy(self):
        return copy.deepcopy(self.positions)

    def radius_copy(self):
        return self.radius.copy()

    def radius_deepcopy(self):
        return copy.deepcopy(self.radius)



    def update_labels(self, tlabel):
        if len(self.time_labels) > 0 :
            self.current_time_index=self.time_labels.index(self.last_read_time_label)
        self.last_read_time_label=list(tlabel)
        if tlabel not in self.time_labels:
            self.time_labels.append(tlabel)
#                    self.time_labels.sort() # normally useless, as it should be already ordered


    def convert_input_pos(self,line):
        values=split(line)
        
        switch=0

        if(len(values)>0):
            if(values[0] == '#'):
                if self.is_file:
                    tlabel=[float(values[1]),self.instream.tell()]
                else:
                    tlabel=[float(values[1]),0]

                self.update_labels(tlabel)
                switch=1
#                print " Read positions ", values[1]

            if len(values) > 10 : # exact number to be fixed. Now it is different in 2d an 3d, due to different output of angular position
                i=int(values[0])
#                self.positions[i]=[float(values[j])+0.5*self.cell_size for j in range(1,4)]
                self.positions[i]=[float(values[j]) for j in range(2,5)]
                self.radius[i]=float(values[1])

        return switch


# end of convert_input(input_line)



    def get_snapshot(self, verbose=False):

        self.old_positions=self.positions.copy()
        self.positions.clear()
        switch=0

        count=0

        while True:
            input_line=self.instream.readline()
            if len(input_line)==0:
                break              #EOF
            if verbose:
                print input_line

            switch+=self.convert_input_pos(input_line)

            if switch==1:
                return True

            
        return False


    def computePairSeparations(self):
        flat_pos = np.ravel(self.positions.values())
        self.pair_sep = np.asarray([ flat_pos - np.roll(flat_pos, 3*i, axis=0) for i in self.range()] )

        # define the shifts to apply if a boundary is crossed
        xcross_shift = np.tile([self.Lx(), 0, 0], self.np()) # x
        ycross_shift = np.tile([0, self.Ly(), 0], self.np()) # y
        xshift = (self.time()-int(self.time()))*self.Lx()    # z (Lees-Edwards)
        zcross_shift = np.tile([xshift, 0, self.Lz()], (self.np(),self.np())) # z

        # easier to start by delta_z

        # case delta_z > 0.5*Lz:
        # determine which elements have to be updated:
        # if delta_z > 0.5*Lz, delta_x and delta_z should be updated, delta_y remains
        # so the 3 entries will appear in the mask as [ True, False, True ]
        zcrossed = np.repeat(self.pair_sep[:,2::3] > 0.5*self.Lz(),2, axis=1) #  [ True, True ]
#        mask = np.insert(zcrossed, np.tile(np.arange(1,2*zcrossed.shape[0],2), (200, 1)), False, axis=1) # insert False in between
        mask = np.insert(zcrossed, np.arange(1,2*zcrossed.shape[0],2), False, axis=1) # insert False in between
        # replace values using the mask
        self.pair_sep = np.where(mask, self.pair_sep-zcross_shift, self.pair_sep)
        
        # case delta_z < -0.5*Lz:
        # same as above
        zcrossed = np.repeat(self.pair_sep[:,2::3] < -0.5*self.Lz(),2, axis=1)
#        mask = np.insert(zcrossed, np.tile(np.arange(1,2*zcrossed.shape[0],2), (200, 1)), False, axis=1)
        mask = np.insert(zcrossed, np.arange(1,2*zcrossed.shape[0],2), False, axis=1)
        # replace values using the mask
        self.pair_sep = np.where(mask, self.pair_sep+zcross_shift, self.pair_sep)
        

        # case delta_x > 0.5*Lx:
        # delta_x alone should be updated, mask [ True, False, False ]
        mask = np.tile(np.asarray(np.zeros(2*self.np()), dtype=bool), (200, 1)) # [ False, False ]
#        mask = np.insert(mask, np.tile(np.arange(0,2*mask.shape[0],2), (200,1)), self.pair_sep[:,0::3] > 0.5*self.Lx(), axis=1) # [ True, False, False ]
        mask = np.insert(mask, np.arange(0,2*mask.shape[0],2), self.pair_sep[:,0::3] > 0.5*self.Lx(), axis=1) # [ True, False, False ]
        # replace values using the mask
        self.pair_sep = np.where(mask, self.pair_sep-xcross_shift, self.pair_sep)

        # case delta_x < -0.5*Lx:
        mask = np.tile(np.asarray(np.zeros(2*self.np()), dtype=bool), (200, 1)) # [ False, False ]
#        mask = np.insert(mask, np.tile(np.arange(0,2*len(mask),2), (200,1)), self.pair_sep[:,0::3] < -0.5*self.Lx(), axis=1) # [ True, False, False ]
        mask = np.insert(mask, np.arange(0,2*len(mask),2), self.pair_sep[:,0::3] < -0.5*self.Lx(), axis=1) # [ True, False, False ]
        # replace values using the mask
        self.pair_sep = np.where(mask, self.pair_sep+xcross_shift, self.pair_sep)


        # case delta_y > 0.5*Ly:
        # delta_y alone should be updated, mask [ False, True, False ]
        mask = np.tile(np.asarray(np.zeros(2*self.np()), dtype=bool), (200, 1)) # [ False, False ]
#        mask = np.insert(mask, np.tile(np.arange(0,2*len(mask),2), (200,1)), self.pair_sep[:,1::3] > 0.5*self.Ly(), axis=1) # [ False, True, False ]
        mask = np.insert(mask, np.arange(0,2*len(mask),2), self.pair_sep[:,1::3] > 0.5*self.Ly(), axis=1) # [ False, True, False ]
        # replace values using the mask
        self.pair_sep = np.where(mask, self.pair_sep-ycross_shift, self.pair_sep)

        # case delta_y < -0.5*Ly:
        # delta_y alone should be updated, mask [ False, True, False ]
        mask = np.tile(np.asarray(np.zeros(2*self.np()), dtype=bool), (200, 1)) # [ False, False ]
#        mask = np.insert(mask, np.tile(np.arange(0,2*len(mask),2), (200,1)), self.pair_sep[:,1::3] < -0.5*self.Ly(), axis=1) # [ False, True, False ]
        mask = np.insert(mask, np.arange(0,2*len(mask),2), self.pair_sep[:,1::3] < -0.5*self.Ly(), axis=1) # [ False, True, False ]
        # replace values using the mask
        self.pair_sep = np.where(mask, self.pair_sep+ycross_shift, self.pair_sep)

        
        self.pair_dist = np.sum(np.reshape(self.pair_sep, (self.np(), self.np(), 3))**2, axis=2)**(1./2)

    def getPairSeparations(self):
        return self.pair_sep

    cdef void cperiodize(self, double *deltax, double *deltay, double *deltaz):
        cdef double xshift
        xshift = (self.time()-int(self.time()))*self.Lx()

        if deltaz[0] > 0.5*self.Lz():
            deltaz[0] = deltaz[0] - self.Lz()
            deltax[0] = deltax[0] - xshift
        else:
            if deltaz[0] < -0.5*self.Lz():
                deltaz[0] = deltaz[0] + self.Lz()
                deltax[0] = deltax[0] + xshift

        if deltay[0] > 0.5*self.Ly():
            deltay[0] = deltay[0] - self.Ly()
        else:
            if deltay[0] < -0.5*self.Ly():
                deltay[0] = deltay[0] + self.Ly()

        if deltax[0] > 0.5*self.Lx():
            deltax[0] = deltax[0] - self.Lx()
        else:
            if deltax[0] < -0.5*self.Lx():
                deltax[0] = deltax[0] + self.Lx()

    cpdef double rad(self,int i):
        return self.radius[i]

    cpdef list pos(self,int i):
        return self.positions[i]

    cpdef list periodized_pos(self,int i):

        cdef double x,y,z
        x=self.positions[i][0]
        y=self.positions[i][1]
        z=self.positions[i][2]
        self.cperiodize(&x, &y, &z)
        return [x, y, z]

    cpdef list pos_diff_vec(self,a,b):
        cdef double deltax
        cdef double deltay
        cdef double deltaz
        cdef double deltatot

        deltax=a[0]-b[0]
        deltay=a[1]-b[1]
        deltaz=a[2]-b[2]
        self.cperiodize(&deltax, &deltay, &deltaz)
        deltatot= libc.math.sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz)
        return [deltax, deltay, deltaz, deltatot]

    cdef void cpos_diff_vec(self, double *a, double *b, double *delta):

        delta[0]=a[0]-b[0]
        delta[1]=a[1]-b[1]
        delta[2]=a[2]-b[2]
        self.cperiodize(&delta[0], &delta[1], &delta[2])
        delta[3]= libc.math.sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2])

    cpdef list pos_diff(self, int i, int j):
        return self.pos_diff_vec(self.positions[i],self.positions[j])

    cdef void cpos_velocity_diff(self, int i, int j, double *deltax, double *deltay, double *deltaz, double *deltatot, double *velox, double *veloy, double *veloz, double *velotot):


        cdef int s
        cdef double old_posi [3]
        cdef double new_posi [3]
        cdef double old_posj [3]
        cdef double new_posj [3]
        for s in range(3):
            old_posi[s]=self.old_positions[i][s]
            new_posi[s]=self.positions[i][s]
            old_posj[s]=self.old_positions[j][s]
            new_posj[s]=self.positions[j][s]

        cdef double old_posdiff [4]
        cdef double new_posdiff [4]
        cdef double dt


        self.cpos_diff_vec(old_posi,old_posj,old_posdiff)
        self.cpos_diff_vec(new_posi,new_posj,new_posdiff)


        dt=self.time()-self.time(-1)

        velox[0]=(new_posdiff[0]-old_posdiff[0])/dt
        veloy[0]=(new_posdiff[1]-old_posdiff[1])/dt
        veloz[0]=(new_posdiff[2]-old_posdiff[2])/dt
        velotot[0]=(new_posdiff[3]-old_posdiff[3])/dt
        deltax[0]=new_posdiff[0]
        deltay[0]=new_posdiff[1]
        deltaz[0]=new_posdiff[2]
        deltatot[0]=new_posdiff[3]


    cpdef list pos_velocity_diff(self, int i, int j):


        cdef int s
        cdef double old_posi [3]
        cdef double new_posi [3]
        cdef double old_posj [3]
        cdef double new_posj [3]
        for s in range(3):
            old_posi[s]=self.old_positions[i][s]
            new_posi[s]=self.positions[i][s]
            old_posj[s]=self.old_positions[j][s]
            new_posj[s]=self.positions[j][s]

        cdef double old_posdiff [4]
        cdef double new_posdiff [4]
        cdef double dt


        self.cpos_diff_vec(old_posi,old_posj,old_posdiff)
        self.cpos_diff_vec(new_posi,new_posj,new_posdiff)
        dt=self.time()-self.time(-1)

        velo=[]
        delta=[]
        for s in range(4):
            velo.append((new_posdiff[s]-old_posdiff[s])/dt)
            delta.append(new_posdiff[s])

        return [delta, velo]


    cpdef list verbose_pos_velocity_diff(self, int i, int j):


        cdef int s
        cdef double old_posi [3]
        cdef double new_posi [3]
        cdef double old_posj [3]
        cdef double new_posj [3]
        for s in range(3):
            old_posi[s]=self.old_positions[i][s]
            new_posi[s]=self.positions[i][s]
            old_posj[s]=self.old_positions[j][s]
            new_posj[s]=self.positions[j][s]

        cdef double old_posdiff [4]
        cdef double new_posdiff [4]
        cdef double dt


        self.cpos_diff_vec(old_posi,old_posj,old_posdiff)
        self.cpos_diff_vec(new_posi,new_posj,new_posdiff)
        dt=self.time()-self.time(-1)

        print "old pos diff " , old_posdiff[0], old_posdiff[1], old_posdiff[2], old_posdiff[3]
        print "new pos diff " , new_posdiff[0], new_posdiff[1], new_posdiff[2], new_posdiff[3]

        velo=[]
        delta=[]
        for s in range(4):
            velo.append((new_posdiff[s]-old_posdiff[s])/dt)
            delta.append(new_posdiff[s])

        print "velo", velo
        return [delta, velo]



    cpdef list instant_relative_velocity(self, int i, int j):

        od=self.pos_diff_vec(self.old_positions[i],self.old_positions[j])
        nd=self.pos_diff_vec(self.positions[i],self.positions[j])
        v=[0.,0.,0.,0.]
        dt=self.time()-self.time(-1)
        for s in range(4):
            v[s]=(nd[s]-od[s])/dt
        return v

    cpdef list average_relative_velocity(self, int i, int j,avg_time, previous_position):
        od=self.pos_diff_vec(previous_position[i],previous_position[j])
        nd=self.pos_diff(i,j)
#        print i, j, od, nd

        v=[0.,0.,0.,0.]
        dt=avg_time
        for s in range(4):
            v[s]=(nd[s]-od[s])/dt
        return v

#end of pos_diff

    cpdef list velocity(self, int i):
        is_computed=True

        for j in [0, -1, -2]:
            if (self.time(j)-libc.math.floor(self.time(j))) < 1e-6:
                    is_computed=False

        if is_computed:
        #        print pos_diff_vec(positions[i],old_positions[i])
            dt=self.time(0)-self.time(-1)
            dr=self.pos_diff_vec(self.positions[i],self.old_positions[i])
            dxinf=self.positions[i][1]*dt # == gamma_dot*y*dt
            dr[0]-=dxinf
            return [ dr[j]/dt for j in range(4) ]
        else:
            return [0.,0.,0.,0.]

    def displacement(self, i):
        dr=self.pos_diff_vec(self.positions[i],self.old_positions[i])
        return dr

    def non_affine_displacement(self, i):
        dr=self.displacement(i)
        dt=self.time(0)-self.time(-1)
        dr[0]-=self.old_positions[i][2]*dt
        dr[3]=libc.math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])
        return dr
        
    def goto_time(self, time):
        if self.is_file:
            for tl in self.time_labels:
                if tl[0]>=time:
                    self.instream.seek(tl[1])
                    self.current_time_index=self.time_labels.index(tl)
                    self.last_read_time_label=tl
                    return
        else:
            sys.stderr.write("Cannot go to time %s".str(time))
            sys.stderr.write("Input is %s".str(instream))
            sys.exit(1)

    def rewind(self): # go to beginning of the file
        if self.is_file:
            self.instream.seek(0)
            self.reset_members()
            return
        else:
            sys.stderr.write("Cannot rewind")
            sys.stderr.write("Input is %s".str(instream))
            sys.exit(1)


    cpdef double time(self, t=0):
        tl_loc=self.current_time_index+t
        return self.time_labels[tl_loc][0]

    cpdef int time_index(self):
        return self.current_time_index


    def part_nb(self):
        return len(self.positions)

    def np(self):
        return self.N

    def rho(self):
        return self.N/self.V

    def vol(self):
        return self.V

    def Lx(self):
        return self.lx

    def Ly(self):
        return self.ly

    def Lz(self):
        return self.lz

    def range(self, int u=0):
        return range(u,self.part_nb())
    #    def range(self, u=1):
    #    return range(u,self.part_nb()+1)

    def is_present(self, i):
        return (i in self.positions)
    
    def dimension(self):
        return self.dim
