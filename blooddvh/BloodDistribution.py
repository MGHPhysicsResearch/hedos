import numpy as np
import pandas as pd
import time

from functools import reduce
from matplotlib.cbook import flatten
from blooddvh import CompartmentModel

"""
Blood distributions over compartment organ
"""
class BloodDistribution:
    __slots__ = ["df", "name", "volume", "dt", "tt", "ttd", "mtt"]
    def __init__(self):
        self.tt   = []
        self.ttd  = []
        self.mtt  = []
        self.dt   = 0
        self.df   = None

    def generate_from_markov(self, markov, names, volume, dt, nb_samples, nb_steps, seed=0):
        """
        Generate a blood distribution from a pre-built markov chain, e.g., from Compartment model
        """
        self.df     = pd.DataFrame(0, index=np.array(range(nb_samples)), columns=np.array(range(nb_steps)))
        self.name   = names
        self.volume = volume
        self.dt     = dt
        
        t  = time.process_time()
        for i in range(nb_samples) :
            start           = self.name[ self.volume.searchsorted( np.random.uniform() ) ]
            self.df.iloc[i] = markov.walk(nb_steps, start, output_indices=True) #seed can be assigned too.
            elapsed_time    = time.process_time() - t
        print("time to generate blood distribution", elapsed_time)

    def save_blood_distribution(self, excel_file, tab_name):
        """
        Save blood path to excel
        - 
        """
        writer = pd.ExcelWriter(excel_file) 
        self.df.to_excel(writer,tab_name)
        master = pd.DataFrame()
        #master.name
        #master.volume
        #mater.dt
        # don't know how to add one row for master tab
        writer.save()
        
    def read_from_excel(self, f_name, s_name, h_name="master"):
        """
        Read blood distribution from Excel file
        TODOs: volume setup
        """
        self.df   = pd.read_excel(f_name, sheet_name=s_name)
        header    = pd.read_excel(f_name, sheet_name=h_name, header=None, index_col=0)
        self.name = list(header.loc["name"])
        self.dt   = header.loc["dt_sec",1]
        #self.volume

    def read_from_df(self, df, name, volume, dt):
        """
        Set up blood distribution from DataFrame
        """
        self.df   = df
        self.name = name
        self.dt   = dt
        self.volume = volume

    def count_consecutives(self, compartment_id):
        """
        Count number of consecutives per particle
        From [1, 1, 2, 2, 1, 1, 1, 1, 9, 1]
        To   [2, 4, 1] -> enters a compartment "1" for 3 times and spent there 2 * dT, 4*dT, dT
        """
        bp_in_compartment = []
        for i, row in self.df.iterrows():
            c = self.df.iloc[i].values == compartment_id 
            np.concatenate(([c[0]], c[:-1] != c[1:], [True]))
            d = np.diff(np.where(np.concatenate(([c[0]], c[:-1] != c[1:], [True])))[0])[::2]
            if np.size(d) > 0 : bp_in_compartment.append(d)
        return bp_in_compartment

    def transition_time(self):
        """
        Calculate transition time (tt), tt-distribution, and mean tt (mtt)
        """
        for i in self.name :
            t = self.count_consecutives(self.name.index(i))
            l = list(flatten(t))
            self.tt.append( t )
            self.ttd.append( l )
            self.mtt.append( [np.mean(l), np.std(l)] )
