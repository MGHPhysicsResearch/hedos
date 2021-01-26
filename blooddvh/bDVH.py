from blooddvh import BloodDistribution
from blooddvh import tDVH
import pandas as pd
import numpy  as np
import copy   as copy

"""
Blood dose calculations using Blood distribution and time-dependent DVH
"""
class bDVH:
    __slots__ = ["dose", "dose_per_fraction", "fractions", "seed", "blood_path", "dT"]
    def __init__(self, blood_path, dT, seed=0):
        """
        blood_path: Expected to be Pandas' Dataframe
        dT        : Time (sec) per step in blood_path
        """
        self.blood_path = blood_path
        self.dT         = dT
        self.seed       = seed

        # total dose
        self.dose          = np.zeros(self.blood_path.shape[0])

        # fraction doses
        # row: BP
        # col: Dose at fraction (it's not accumulation dose)
        self.dose_per_fraction = np.zeros((self.blood_path.shape[0],1))

        # fraction size
        self.fractions = []

    def volume_gt_dose(self, threshold):
        """
        threshold: threshold dose (Gy)
        return volume (%) of BP g.t threshold, mean dose, std dose
        """
        # return mean/std of dose?
        bp     = self.dose[self.dose  >  threshold]
        bp_vol = 100.0 * bp.size/self.dose.size
        return bp_vol, bp.mean(), bp.std()

    def dose_at_top_volume(self, threshold):
        """
        threshold (%) : upper volume, e.g., 2 for 2 %
        return dose (Gy) at volume, mean dose, std dose
        """
        np_top = int(self.blood_path.shape[0] * threshold * 0.01)
        t = np.sort(self.dose)
        return t[-np_top], t[-np_top:].mean(), t[-np_top:].std()
        
    # Adding compartment DVHs per courses
    # For multiple courses plan, i.e., different fraction sizes, this methods should get called multiple times.
    def add_dose(self, beams, compartment_ids, start_time = 0.0, fraction = 1):
        """
        beams: a list  or one tDVH 
        compartment_ids: a list or one compartment id. number should match to beams
        fractions  : number of fraction. This will expand the column of self.dose 
        start_time : time step when 
        """
        
        dose_1st = np.zeros(self.blood_path.shape[0])

        # First fraction dose
        if isinstance(beams, list) and isinstance(compartment_ids, list) :
            # in case, multiple-beams and -compartments per course
            for i, x in enumerate(zip(beams, compartment_ids)):
                dose_1st += self.__add_dose(x[0], x[1], start_time)
        else:
            # only one beam & compartment is given
            dose_1st = self.__add_dose(beams, compartment_ids, start_time)

        # Let's repeat 1st-dose N times
        # and update total dose
        self.__repeat(dose_1st, fraction)

    # shuffle and add 1st-fraction dose by fx times
    def __repeat(self, dose_1st, fractions = 1 ):
        """
        dose_1st: first dose (row) for blood particles
        fx      : number of fractions 
        """

        # Add 1st dose to dose_per_fraction data
        if np.size(self.fractions) == 0 :
            self.dose_per_fraction[0:,0] = dose_1st
        else:
            self.dose_per_fraction = np.hstack((self.dose_per_fraction, np.transpose([dose_1st])))

        # Fraction size update
        self.fractions.append(fractions)

        # Shuffle & Repeat
        if fractions > 1 :
            # copy 1st dose (N -1 ) times because 1st dose is already appended above
            d_group = np.vstack([dose_1st] * (fractions-1) ).transpose()
            for i in range(0, fractions-1):
                # shuffle to make fraction dose
                np.random.shuffle(d_group[0:,i])
                
            # after for-loop, update dose_per_fraction table
            self.dose_per_fraction = np.hstack((self.dose_per_fraction, d_group))

        # Update final dose
        self.dose += np.sum(self.dose_per_fraction[0:, -fractions:], axis=1)
        
        
    # Beam may have multiple DVHs
    # One compartment model
    # Allowed from internal call
    def __add_dose(self, beam, compartment_id, start_time=0.0):
        """
        beam::tDVH : time-dependent DVHs
        compartment_id : compartment gets dose
        beam_start : time-point beam starts to delivery
        """
        d = np.zeros(self.blood_path.shape[0])
        assert( start_time + beam.dvh[0][-1] <= self.blood_path.shape[1] * self.dT )
        
        # For index, time, and DVH function
        for i, (_, dvh) in enumerate(zip(beam.dvh[0], beam.dvh[1])):
            if dvh == None:
                continue
        
            # Two time-points for blood path (t0,t1) get irradiated by dose.
            t0   = start_time if i == 0 else start_time + beam.dvh[0][i-1]
            t1   = start_time + beam.dvh[0][i]
            dt1_t0 = t1 - t0

            # Calculate index of blood-particle distribution (BPD)
            lower_bound = np.floor(t0/self.dT)
            upper_bound = np.ceil(t1/self.dT)
            idx0        = int(lower_bound)
            idx1        = int(upper_bound)
            distance    = idx1 - idx0
            #print("i, t0, t1, idx0, idx1, distance: ", t0, t1, idx0, idx1, distance)
    
            for j, r in enumerate(self.blood_path.iloc[0:, idx0:idx1].values):
        
                # 1st column: add dose before loop to account partial delivery
                if r[0] == compartment_id :
                    d[j] += (self.dT - (t0 - lower_bound)) * dvh(np.random.uniform())/dt1_t0

                # path from second to last
                for _ in range(np.count_nonzero(r[1:] == compartment_id)):
                    d[j] += self.dT*dvh(np.random.uniform())/dt1_t0

                # last column: minus dose in case partial delivery
                if r[distance-1] == compartment_id :
                    d[j] -= (upper_bound - t1) * dvh(np.random.uniform())/dt1_t0

        # outside for-loop
        return d
