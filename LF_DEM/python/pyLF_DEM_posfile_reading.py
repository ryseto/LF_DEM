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

        if self.ly == 0.: # 2d case
            self.V = self.lx*self.lz
            self.dim=2
        else:               # 3d case
            self.V = self.lx*self.ly*self.lz
            self.dim=3

        self.reset_members()

    def reset_members(self):
        self.positions = np.zeros((self.N, 3))
        self.old_positions = np.zeros((self.N, 3))
        self.radius = np.zeros(self.N)
        self.last_read_time_label=[]
        self.time_labels=[]
        self.current_time_index=0
        test_input=self.get_snapshot()
        self.first_time=self.time()
        if test_input==False:
            print "Invalid input "
            sys.exit(1)
        

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
                self.positions[i,:] = np.array([values[j] for j in range(2,5)], dtype=float)
                self.radius[i] = float(values[1])

        return switch


    def get_snapshot(self, verbose=False):

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

    def periodize(self,delta):
        
        deltax = delta[0]
        deltay = delta[1]
        deltaz = delta[2]
        xshift = (self.time()-int(self.time()))*self.lx

        if deltaz > 0.5*self.lz:
            deltaz = deltaz - self.lz
            deltax = deltax - xshift
        else:
            if deltaz < -0.5*self.lz:
                deltaz = deltaz + self.lz
                deltax = deltax + xshift
                
        if deltay > 0.5*self.ly:
            deltay = deltay - self.ly
        else:
            if deltay < -0.5*self.ly:
                deltay = deltay + self.ly

        if deltax > 0.5*self.lx:
            deltax = deltax - self.lx
        else:
            if deltax < -0.5*self.lx:
                deltax = deltax + self.lx

        return [deltax, deltay, deltaz]


    def computePairSeparations(self):
        self.pair_sep = self.positions[:,None,:] - self.positions

        # define the shifts to apply if a boundary is crossed
        xshift = (self.time()-int(self.time()))*self.lx    # z (Lees-Edwards)

        # easier to start by delta_z

        # case delta_z > 0.5*Lz:
        # determine which elements have to be updated:
        crossed = self.pair_sep[:,:,2] > 0.5*self.lz
        # shift values
        self.pair_sep = np.where(crossed[:,:,None], self.pair_sep-[xshift, 0, self.lz], self.pair_sep)
        
        # case delta_z < -0.5*Lz:
        # same as above
        crossed = self.pair_sep[:,:,2] < -0.5*self.lz
        # shift values
        self.pair_sep = np.where(crossed[:,:,None], self.pair_sep+[xshift, 0, self.lz], self.pair_sep)

        # case delta_x > 0.5*Lx:
        crossed = self.pair_sep[:,:,0] > 0.5*self.lx
        while crossed.any():
            self.pair_sep = np.where(crossed[:,:,None], self.pair_sep-[self.lx, 0, 0], self.pair_sep)
            crossed = self.pair_sep[:,:,0] > 0.5*self.lx

        # case delta_x < -0.5*Lx:
        crossed = self.pair_sep[:,:,0] < -0.5*self.lx
        while crossed.any():
            self.pair_sep = np.where(crossed[:,:,None], self.pair_sep+[self.lx, 0, 0], self.pair_sep)
            crossed = self.pair_sep[:,:,0] < -0.5*self.lx

        # case delta_y > 0.5*Ly:
        crossed = self.pair_sep[:,:,1] > 0.5*self.ly
        self.pair_sep = np.where(crossed[:,:,None], self.pair_sep-[0, self.ly, 0], self.pair_sep)

        # case delta_y < -0.5*Ly:
        crossed = self.pair_sep[:,:,1] < -0.5*self.ly
        self.pair_sep = np.where(crossed[:,:,None], self.pair_sep+[0, self.ly, 0], self.pair_sep)

        self.pair_dist = np.sqrt((self.pair_sep**2).sum(axis=2))


    def periodized_pos(self,i):
        return self.periodize(self.positions[i])



    def displacement(self, i):
        dr=self.pos_diff_vec(self.positions[i],self.old_positions[i])
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


    def time(self, t=0):
        tl_loc=self.current_time_index+t
        return self.time_labels[tl_loc][0]

    def rho(self):
        return self.N/self.V

    def range(self, u=0):
        return range(u,self.N)
