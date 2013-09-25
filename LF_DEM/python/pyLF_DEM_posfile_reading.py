import math
from string import *
import sys
import copy
import numpy as np

class Pos_Stream:

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

    def reset_members(self):
        self.positions_dict=dict()
        self.old_positions=dict()
        self.radius_dict=dict()
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
                self.positions_dict[i]=[float(values[j]) for j in range(2,5)]
                self.radius_dict[i]=float(values[1])

        return switch


    def get_snapshot(self, verbose=False):

#        self.old_positions=self.positions.copy()
        self.positions_dict.clear()
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
                self.positions = np.asarray(self.positions_dict.values())
                return True

            
        return False

    def periodize(self,delta):
        
        deltax = delta[0]
        deltay = delta[1]
        deltaz = delta[2]
        xshift = (self.time()-int(self.time()))*self.Lx()

        if deltaz > 0.5*self.Lz():
            deltaz = deltaz - self.Lz()
            deltax = deltax - xshift
        else:
            if deltaz < -0.5*self.Lz():
                deltaz = deltaz + self.Lz()
                deltax = deltax + xshift
                
        if deltay > 0.5*self.Ly():
            deltay = deltay - self.Ly()
        else:
            if deltay < -0.5*self.Ly():
                deltay = deltay + self.Ly()

        if deltax > 0.5*self.Lx():
            deltax = deltax - self.Lx()
        else:
            if deltax < -0.5*self.Lx():
                deltax = deltax + self.Lx()

        return [deltax, deltay, deltaz]

    def pos(self,i):
        i=str(i)
        return self.positions_dict[i]

    def rad(self,i):
        i=str(i)
        return self.radius_dict[i]



    def computePairSeparations(self):
        flat_pos = np.ravel(self.positions)
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

        self.pair_sep = np.reshape(self.pair_sep, (self.np(), self.np(), 3))
        
        self.pair_dist = np.sum(self.pair_sep**2, axis=2)**(1./2)


    def periodized_pos(self,i):
        i=str(i)
        return self.periodize(self.positions_dict[i])


    def pos_diff_vec(self,a,b):

        
        deltax=a[0]-b[0]
        deltay=a[1]-b[1]
        deltaz=a[2]-b[2]
        print "pos ", a
        print "pos ", b
        print [deltax, deltay, deltaz]
        [deltax, deltay, deltaz]=self.periodize([deltax, deltay, deltaz])
        print [deltax, deltay, deltaz]
        deltatot= math.sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz)
        return [deltax, deltay, deltaz, deltatot]
    
    def pos_diff(self,i,j):
        i=str(i)
        j=str(j)

        return self.pos_diff_vec(self.positions_dict[i],self.positions_dict[j])


    def displacement(self, i):
        dr=self.pos_diff_vec(self.positions_dict[i],self.old_positions_dict[i])
        return dr

    def non_affine_displacement(self, i):
        dr=self.displacement(i)
        dt=self.time(0)-self.time(-1)
        dr[0]-=self.old_positions_dict[i][1]*dt
        dr[3]=math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])
        return dr

    def instant_relative_velocity(self,i,j):
        i=str(i)
        j=str(j)

        od=self.pos_diff_vec(self.old_positions_dict[i],self.old_positions_dict[j])
        nd=self.pos_diff_vec(self.positions_dict[i],self.positions_dict[j])
        v=[0.,0.,0.,0.]
        dt=self.time()-self.time(-1)
        for s in range(4):
            v[s]=(nd[s]-od[s])/dt
        return v

    def average_relative_velocity(self,i,j,avg_time, previous_position):
        i=str(i)
        j=str(j)
        
        od=self.pos_diff_vec(previous_position[i],previous_position[j])
        nd=self.pos_diff(i,j)
#        print i, j, od, nd

        v=[0.,0.,0.,0.]
        dt=avg_time
        for s in range(4):
            v[s]=(nd[s]-od[s])/dt
        return v
            
#end of pos_diff

    def velocity(self,i):
        
        i=str(i)
        #        print pos_diff_vec(positions_dict[i],old_positions_dict[i]) 
        dt=self.time(0)-self.time(-1)
        dr=self.pos_diff_vec(self.positions_dict[i],self.old_positions_dict[i])
        dxinf=self.positions_dict[i][1]*dt # == gamma_dot*y*dt
        dr[0]-=dxinf
        return [ dr[j]/dt for j in range(4) ]


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


    def time(self, t=0):
        tl_loc=self.current_time_index+t
        return self.time_labels[tl_loc][0]


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

    def range(self, u=0):
        return range(u,self.np())
#     #    def range(self, u=1):
#     #    return range(u,self.part_nb()+1)

    # def is_present(self, i):
    #     return (i in self.positions_dict)
    
    def dimension(self):
        return self.dim
