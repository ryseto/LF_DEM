import math
from string import *
import sys
import copy

class Pos_Stream:
    
    
    def __init__(self, input_stream, N, phi):

        self.instream=input_stream
        try:
            self.instream.tell()
            self.is_file=True
        except IOError:
            self.is_file=False

        V=4.*math.acos(-1)*N/(3.*phi)
        self.cell_size=math.pow(V, 1./3.)

        self.reset_members()
        
    def __iter__(self):
        return self.positions.__iter__()

    def reset_members(self):
        self.positions=dict()
        self.radius=dict()
        self.GUh_stress=dict()
        self.GUc_stress=dict()
        self.xFc_stress=dict()
        self.old_positions=dict()
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
                
            if len(values) > 10 : # exact number to be fixed. Now it is different in 2d an 3d, due to different output of angular position
                i=values[0]
#                self.positions[i]=[float(values[j])+0.5*self.cell_size for j in range(1,4)]
                self.positions[i]=[float(values[j]) for j in range(2,5)]
                self.radius[i]=float(values[1])
                self.GUh_stress[i]=float(values[11])
                self.xFc_stress[i]=float(values[12])
                self.GUc_stress[i]=float(values[13])
        return switch
            
# end of convert_input(input_line)



    def get_snapshot(self, verbose=False):

        self.old_positions=self.positions.copy()
        self.positions.clear()
        switch=0
        
        count=0
        input_line='start'
        while input_line!='':

            input_line=self.instream.readline()
            if verbose:
                print input_line

            switch+=self.convert_input_pos(input_line)

            if switch==1:
                return True

        return False

    def periodize(self,delta):

        xshift = self.time()
        deltax=delta[0]
        deltay=delta[1]
        deltaz=delta[2]


        if deltaz > 0.5*self.cell_size:
            deltaz = deltaz - self.cell_size
            deltax = deltax - xshift
        else:
            if deltaz < -0.5*self.cell_size:
                deltaz = deltaz + self.cell_size
                deltax = deltax + xshift

        if deltay > 0.5*self.cell_size:
            deltay = deltay - self.cell_size
        else:
            if deltay < -0.5*self.cell_size:
                deltay = deltay + self.cell_size

        if deltax > 0.5*self.cell_size:
            deltax = deltax - self.cell_size
        else:
            if deltax < -0.5*self.cell_size:
                deltax = deltax + self.cell_size

        return [deltax, deltay, deltaz]

    def pos(self,i):
        i=str(i)
        return self.positions[i]

    def rad(self,i):
        i=str(i)
        return self.radius[i]

    def periodized_pos(self,i):
        i=str(i)
        return self.periodize(self.positions[i])


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

        return self.pos_diff_vec(self.positions[i],self.positions[j])


    def displacement(self, i):
        dr=self.pos_diff_vec(self.positions[i],self.old_positions[i])
        return dr

    def non_affine_displacement(self, i):
        dr=self.displacement(i)
        dt=self.time(0)-self.time(-1)
        dr[0]-=self.old_positions[i][1]*dt
        dr[3]=math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])
        return dr

    def instant_relative_velocity(self,i,j):
        i=str(i)
        j=str(j)

        od=self.pos_diff_vec(self.old_positions[i],self.old_positions[j])
        nd=self.pos_diff_vec(self.positions[i],self.positions[j])
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
        #        print pos_diff_vec(positions[i],old_positions[i]) 
        dt=self.time(0)-self.time(-1)
        dr=self.pos_diff_vec(self.positions[i],self.old_positions[i])
        dxinf=self.positions[i][1]*dt # == gamma_dot*y*dt
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

    def part_nb(self):
        return len(self.positions)
    
    def rho(self):
        return float(self.part_nb)/float(pow(self.cell_size,3))

    def cell_lin_size(self):
        return self.cell_size
    def range(self, u=0):
        return range(u,self.part_nb())

    def is_present(self, i):
        return (i in self.positions)

    
