import numpy as np
import pandas as pd
import time

from os import path
from functools import reduce
from matplotlib.cbook import flatten
from blooddvh import CompartmentModel

"""
Blood distributions over compartment organ
has path characterizations, e.g., mean transition time, 
"""
class BloodDistribution:
    __slots__ = ["df", "names", "volumes", "dt", "tv","tt", "ttd", "mtt","rt", "rtd", "mrt", "scales", "shapes"]
    def __init__(self):
        self.tt   = {}
        self.ttd  = {}
        self.mtt  = {}
        self.rt   = {}
        self.rtd  = {}
        self.mrt  = {}
        self.tv   = [] 
        self.dt   = 0
        self.df   = None

    def generate_from_markov(self, markov, names, volumes, dt, nb_samples, nb_steps, seed=0):
        """
        Generate a blood distribution from a pre-built markov chain, e.g., from Compartment model
        """
        self.df      = pd.DataFrame(0, index=np.array(range(nb_samples)), columns=np.array(range(nb_steps)), dtype=np.uint8)
        self.names   = names
        self.volumes = volumes
        self.dt      = dt
        self.scales  = np.zeros(len(self.names)) 
        self.shapes  = np.zeros(len(self.names)) 
        
        t  = time.process_time()
        for i in range(nb_samples) :
            start           = self.names[ self.volumes.searchsorted( np.random.uniform() ) ]
            self.df.iloc[i] = markov.walk(nb_steps, start, output_indices=True) #seed can be assigned too.
        print("time to generate blood distribution", time.process_time()-t)


    def generate_from_markov_weibull(self, wbmc, names, volumes, dt, nb_samples, nb_steps, seed=0):
        """
        Generate a blood distribution from a pre-built markov chain using Weibull distribution, e.g., from Compartment model
        """
        self.df      = pd.DataFrame(0, index=np.array(range(nb_samples)), columns=np.array(range(nb_steps)), dtype=np.uint8)
        self.names   = names
        self.volumes = volumes
        self.dt      = dt
        self.scales  = np.zeros(len(self.names)) 
        self.shapes  = np.zeros(len(self.names)) 
        for i,v in enumerate(wbmc.comp):
            self.scales[i] = v.scale
            self.shapes[i] = v.shape
        
        t  = time.process_time()
        for i in range(nb_samples) :
            compartment_id  = self.volumes.searchsorted( np.random.uniform() )
            compartment_t0  = wbmc.comp[compartment_id].cdf_inv( np.random.uniform() )
            self.df.iloc[i,0]  = compartment_id 
            self.df.iloc[i,1:] = wbmc.walk(nb_steps-1, compartment_id, compartment_t0, self.dt)

        print("time to generate blood distribution", time.process_time()-t)

    def save(self, excel_file, tab_name, ext='bin'):
        """
        Save blood path to excel & tab name
        - 
        """
        fmode   = 'a' if path.exists(excel_file) else 'w'
        writer = pd.ExcelWriter(excel_file, mode=fmode)
        
        master = pd.DataFrame(index=['names', 'volumes','scales','shapes', 'dt_sec'], columns=range(len(self.volumes)))
        master.loc['names']        = self.names
        master.loc['volumes']      = self.volumes
        master.loc['scales']       = self.scales
        master.loc['shapes']       = self.shapes
        master.loc['dt_sec', 0]    = self.dt
        master.loc['samples', 0]   = self.df.shape[0]
        master.loc['steps', 0]     = self.df.shape[1]
        master.loc['path_file', 0] = tab_name + "." + ext
        i = -1*excel_file[::-1].find("/")
        self.df.values.astype('uint8').tofile( excel_file[:i] + tab_name + "." + ext )
        
        master.to_excel(writer, tab_name)
        
        writer.save()
        
    def read_from_excel(self, f_name, t_name):
        """
        Read blood distribution from Excel file
        
        """
        header       = pd.read_excel(f_name, sheet_name=t_name, index_col=0, engine="openpyxl")
        self.names   = list(header.loc["names"])
        self.volumes = list(header.loc["volumes"])
        self.scales  = list(header.loc["scales"])
        self.shapes  = list(header.loc["shapes"])
        self.dt      = header.loc["dt_sec",0]
        i = -1*f_name[::-1].find("/")
        bp_file = f_name[:i] + header.loc["path_file",0]
        print(bp_file)
        bp_raw  = np.fromfile(bp_file, dtype=np.uint8).reshape([header.loc["samples",0], header.loc["steps",0]])
        self.df = pd.DataFrame(bp_raw)
        
    def count_consecutives(self, compartment_id):
        """
        Count number of consecutives per particle from https://stackoverflow.com/a/24343375
        From [1, 1, 2, 2, 1, 1, 1, 1, 9, 1]
        To   [2, 4, 1] -> enters a compartment "1" for 3 times and spent there 2 * dT, 4*dT, dT
        inverse : to calculate distance between enters 
        """
        bp_in_compartment = []
        
        for i, row in self.df.iterrows():
            
            c = self.df.iloc[i].values == compartment_id 
            np.concatenate(([c[0]], c[:-1] != c[1:], [True]))
            d = np.diff(np.where(np.concatenate(([c[0]], c[:-1] != c[1:], [True])))[0])[::2]
            
            if np.size(d) > 0 : bp_in_compartment.append(d)
        return bp_in_compartment

    def count_distances(self, compartment_id):
        """
        Count distances between a compartment.
        From [1, 1, 2, 2, 1, 1, 1, 1, 2, 9, 1]
        To   [2, 2] ->  a BP returns to compartment "1" for 2 times with their distance of 2, respectively.
        """
        distances_compartment = []

        for i, row in self.df.iterrows():
            c = row.values != compartment_id 
            np.concatenate(([c[0]], c[:-1] != c[1:], [True]))
            d = np.diff(np.where(np.concatenate(([c[0]], c[:-1] != c[1:], [True])))[0])[::2]
            
            if row.values[0]  != compartment_id :
                d = d[1:]
            if row.values[-1] != compartment_id :
                d = d[:-1]
            if np.size(d) > 0 : distances_compartment.append(d)
            
        return distances_compartment

    def transition_time(self, names=[]):
        """
        Calculate transition time (tt), tt-distribution(ttd), and mean tt (mtt)
        """
        if len(names) == 0 :
            names = self.names
        
        for i in names :
            t = self.count_consecutives(self.names.index(i))
            l = list(flatten(t))
            self.tt[i]  = t
            self.ttd[i] = l
            self.mtt[i] = [np.mean(l), np.std(l)]

    def recurrence_time(self, names=[]):
        """
        Calculate recurrence time (rt), rt-distribution(rtd), and mean rt (mrt)
        """
        if len(names) == 0 :
            names = self.names

        for i in names :
            t = self.count_distances(self.names.index(i))
            l = list(flatten(t))
            self.rt[i]  = t
            self.rtd[i] = l
            self.mrt[i] = [np.mean(l), np.std(l)]
            
    def temporal_volume(self):
        """
        Calculate volume changes in time. # of BP x # of time-steps
        """
        self.tv = np.zeros([len(self.names), self.df.shape[1]])

        for t in range(self.df.shape[1]):
            # loop over time-steps
            for (c,v) in self.df[t].value_counts(sort=False, normalize=True).items():
                self.tv[c, t] = v

