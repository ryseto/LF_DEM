#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#coding=utf-8

import sys
import numpy
cimport numpy
import string

cdef class TwoTimeCorrelator:

    cdef:
        k
        correlator
        init_positions
        int max_t, init_t


    def __init__(self, wave_vector, positions, init_time, max_time):

        self.k=numpy.array(wave_vector)
        self.init_positions=dict()
        for key in positions:
            self.init_positions[key]=numpy.array(positions[key])
        self.max_t=max_time
        self.init_t=init_time
        self.correlator=[]
        self.update(self.init_positions, init_time)
        
    cpdef update(self, positions, time):
        cdef int key=0
        cdef double cor=0.
        cdef int length=len(positions)
        
        while key < length:
            deltar=numpy.array(positions[key])-self.init_positions[key]
            cor+=numpy.cos(numpy.dot(self.k, deltar))
            key+=1

#        print [time-self.init_t, cor]
        cor/=length
        self.correlator.append([time-self.init_t, cor])


    def alive(self,time):
        if time-self.init_t < self.max_t:
            return True
        else:
            return False
        
    def last_point(self):
        return self.correlator[-1]

    def print_to(self, stream):
        headline='time correlation'+'\n'
        stream.write(headline)

        for point in self.correlator:
            out_string=" ".join(map(str, point))+'\n'
            stream.write(out_string)

    
    

    
class TwoTimeCorrelatorSet:

    def __init__(self, wave_vector, starting_period, max_time):

        self.k=numpy.array(wave_vector)
        self.max_t=max_time
        self.correlators=[]
        self.last_start=-2*starting_period
        self.start_period=starting_period

        self.avg_correlator=dict()
        # i=0
        # while i < max_time:
        #     self.avg_correlator.append([i, []])
        #     i+=



    def add_to_avg(self, time, correl_value):
        if time not in self.avg_correlator:
            self.avg_correlator[time] = []

        self.avg_correlator[time].append(correl_value)

#        self.avg_correlator[bin_lab][1]+=correl_value
#        self.avg_correlator[bin_lab][2]+=1.

        
    def update(self, positions, time):
        
        if (time-self.last_start) > self.start_period:
            new_correlator=TwoTimeCorrelator(self.k, positions, time, self.max_t)
            self.correlators.append(new_correlator)
            self.last_start=time

        if not self.correlators[0].alive(time):
            self.correlators.pop(0)

        for correl in self.correlators:
            correl.update(positions, time)
            lp=correl.last_point()
            self.add_to_avg(lp[0], lp[1])

    def wave_vector(self):
        return self.k

    def get_avg_and_stddev(self, points):
        avg=0.
        stddev=0.

        if len(points)>0:
            for value in points:
                avg+=value 
            avg/=len(points)

            for value in points:
                stddev+=(value-avg)**2
            stddev/=len(points)
            stddev=numpy.sqrt(stddev)
        return [ avg, stddev ]

    def print_to(self, stream):
        headline='k = '+" ".join(map(str, self.k))+'\n'+'time correlation'+'\n'
        stream.write(headline)

        for point in sorted(self.avg_correlator.iterkeys()):
            [avg, stddev]=self.get_avg_and_stddev(self.avg_correlator[point])
            
            out_string=str(point)+' '+str(avg)+' '+str(stddev)+'\n'
#            out_string=str(point[0])+' '+" ".join(map(str, point[1]))+'\n'
#            sys.stdout.write(out_string)
            stream.write(out_string)


    def export_correlator(self):

        cor_to_export = dict()
        for point in sorted(self.avg_correlator.iterkeys()):
            [avg, stddev]=self.get_avg_and_stddev(self.avg_correlator[point])
            cor_to_export[point]=[avg, stddev]

        return cor_to_export
