import numpy as np
import pandas as pd
import time
#from operator import itemgetter
from functools import reduce
from matplotlib.cbook import flatten

class DVHt:
    __slots__ = ["dvh", "seed", "start_time"]

    def __init__(self, seed=0):
        """ 
        Initialize members, dvh, seed, and start_time

        dvh[0] : time, e.g., elapsed time (array) from start
        dvh[1] : dvh function (list)
        """
        self.dvh  = [np.array([]),[]]
        self.seed = seed
        self.start_time = 0.0

    def add(self, beam_on_time, dvh_function):
        """
        Fill dvh functions as a function of time

        Parameters
        beam_on_time : delivery time of this dvh time. Time is relative to beam start
        dvh_function : DVH function (lambda or intpl1d)
        """
        self.dvh[0] = np.append(self.dvh[0], self.dvh[0][-1] + beam_on_time)
        self.dvh[0].append(dvh_function)
    
    def add_csv(self, beam_on_time, dvh_file):
        """
        Parameters
        beam_on_time : delivery time of this dvh time
        dvh_file     : a CSV file with DVH curve
        """
        pass
    
    def set_start_time(self, start_time):
        
        self.start_time = start_time
    
    def get_dose(self, time_from, time_to):
        """
        Return dose between [t0, t1]

        Parameters
        t0, t1: time interval
        """
        # Find dvh
        
        idx = self.dvh[0].searchsorted([time_from, time_to])
        #t   = self.dvh[0][i0], self.dvh[0][i1] #

        for i in idx:
            Ti    = [self.start_time, self.dvh[0][idx[0]] ]
            Ti[0] = self.start_time if idx[0] == 0 else ifself.dvh[0][idx[0]-1] 
        
        # get index of t0
        # get index of t1
        # distance 
        
        v = uniform()
        # d * u : where u: 0-1 for dT
        
        pass
        
