from blooddvh import BloodDistribution
from blooddvh import tDVH
import pandas as pd
import numpy  as np

"""
Blood dose calculations using Blood distribution and time-dependent DVH
"""
class bDVH:
    __slots__ = ["dose", "seed"]
    def __init__(self, seed=0):
        self.seed = seed


    def compute(self, blood_distribution, time_dependent_dvh, organ_id, beam_start=0.0, beam_stop=None):
        """
        Parameter
        blood_distribution::BloodDistribution
        time_dependent_dvh::tDVH
        organ id: dose is accumulated only when a blood particle is inside an organ
        """
        
        # update beam_stop time if it's not given
        beam_stop = beam_start + time_dependent_dvh.dvh[0][-1] if beam_stop == None else beam_stop

        # check beam_stop is outside of blood distribution's time, then produce error 
        assert( beam_stop < (blood_distribution.df.columns[-1] + 1) * blood_distribution.dt )

        # For output
        self.dose = np.zeros( blood_distribution.df.shape[0])

        # [t0, t1] : time window to select blood particles.
        T0 = T1 = beam_start
        for dvh_id in range( np.shape(time_dependent_dvh.dvh)[1]) :
            acc_t  = time_dependent_dvh.dvh[0][dvh_id] #accumulated time
            dvh_f  = time_dependent_dvh.dvh[1][dvh_id] #dvh function in time-interval

            # if no dvh is found in this time-window, skip
            if dvh_f == None :
                continue

            # update upper limit of time-window
            T1 = beam_start + acc_t
            dT = (T1 - T0)

            # print("T0, T1, dT", T0, T1, dT)
            # find column index of blood distribution for given time-window
            c0 = blood_distribution.df.columns.searchsorted(T0, side='right') - 1 
            c1 = blood_distribution.df.columns.searchsorted(T1, side='left') 

            # iterates a subset of DataFrame within time-window
            for i, row in blood_distribution.df.loc[0:, c0:c1].iterrows():

                # map if a blood particle is in organ otherwise false
                is_in_organ = row.values == organ_id

                # number of 'true' but this may overestimate because beam may start in-between blood-step
                N = np.float64(np.sum(is_in_organ))
                #print(N)
                                
                # correction for partial dose if a blood is in organ at first or last index
                if is_in_organ[0] == True :
                    N -= (T0 - c0 * blood_distribution.dt) / blood_distribution.dt 

                if is_in_organ[-1] == True :
                    N -= (c1 * blood_distribution.dt - T1) / blood_distribution.dt 
                
                # call dvh N-times, replace 1.0 with blood-step dt
                for k in range(np.int32(np.floor(N))):
                    self.dose[i] += blood_distribution.dt * (dvh_f(np.random.uniform()) / dT)

                # apply fraction (minus)
                self.dose[i] -= (N - np.floor(N)) * blood_distribution.dt * (dvh_f(np.random.uniform())/dT)
                    
                # update lower limit of time-window for next DVH
                T0 = T1 
