#coding=utf-8

cdef class Histograms:
    cdef:
        histogram
        events

    cdef kadd(self, key, value)
    cdef kappend(self, key, value)
    cdef kreplace(self, key, value)
    

    
