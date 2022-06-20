import numpy as np


"""
Blood dose calculations using Blood distribution and time-dependent DVH
"""
class bDVH:
    __slots__ = ['dose', 'dose_per_fraction', 'fractions', 'seed',
                 'blood_path', 'dT']
    def __init__(self, blood_path, dT, seed=0):
        """

        Parameters
        ----------
        blood_path : pd.DataFrame
            The motion matrix for the simulated blood particles.
        dT : float/int
            The time (seconds) per step in `blood_path`.

        Returns
        -------
        N/A

        """
        self.blood_path = blood_path
        self.dT = dT
        self.seed = seed
        # The total dose
        self.dose = np.zeros(self.blood_path.shape[0])
        # The fraction doses
        # row: number of blood particles
        # col: the dose at the desired fraction (it's not accumulated dose)
        self.dose_per_fraction = np.zeros((self.blood_path.shape[0], 1))
        # The fraction size
        self.fractions = []

    def volume_gt_dose(self, threshold):
        """
        Get the volume of blood particles greater than a specified dose.

        Parameters
        ----------
        threshold : float/int
            The threshold dose (Gy).

        Returns
        -------
        N/A : (float/int, float/int, float/int)
            The volume (%) of blood particles greater than the threshold dose,
            mean dose, and standard deviation dose.

        """
        bp = self.dose[self.dose > threshold]
        bp_vol = 100.0 * bp.size/self.dose.size
        return bp_vol, bp.mean(), bp.std()

    def dose_at_top_volume(self, threshold):
        """
        Get the dose of blood particles greater than a specified volume.

        Parameters
        ----------
        threshold : float/int
            The threshold volume (%).

        Returns
        -------
        N/A : (float/int, float/int, float/int)
            The dose (Gy) of blood particles greater than the threshold
            volume, mean dose, and standard deviation dose.

        """
        np_top = int(self.blood_path.shape[0] * threshold * 0.01)
        t = np.sort(self.dose)
        return t[-np_top], t[-np_top:].mean(), t[-np_top:].std()

    def add_dose(self, beams, compartment_ids, start_time=0.0, fraction=1):
        """
        Add the compartment DVHs per courses. For multiple courses plan, i.e.,
        different fraction sizes, this methods should get called multiple
        times.

        Parameters
        ----------
        beams : tDVH, list[tDVH]
            The time-dependent DVH for the applied beams.
        compartment_ids : int, list[int]
            The compartment ids to add the dose.
        start_time : float/int, optional
            The start time for the applied dose.
        fractions : int, optional
            The number of fractions. This will expand the columns of
            `self.dose`.

        Returns
        -------
        N/A

        """
        dose_1st = np.zeros(self.blood_path.shape[0])
        # The first fraction dose
        if isinstance(beams, list) and isinstance(compartment_ids, list):
            # In case, multiple-beams and compartments per course
            for _,x in enumerate(zip(beams, compartment_ids)):
                dose_1st += self.__add_dose(x[0], x[1], start_time)
        else:
            # Only one beam and compartment is given
            dose_1st = self.__add_dose(beams, compartment_ids, start_time)
        # Repeat 1st-dose N times and update the total dose
        self.__repeat(dose_1st, fraction)

    def __repeat(self, dose_1st, fractions=1):
        """
        Shuffle and add the 1st-fraction dose by fx times.

        Parameters
        ----------
        dose_1st : float/int
            The first dose (row) for the blood particles.
        fractions : int, optional
            The number of fractions.

        Returns
        -------
        N/A

        """
        # Add the 1st dose to `dose_per_fraction` data
        if np.size(self.fractions) == 0:
            self.dose_per_fraction[0:,0] = dose_1st
        else:
            self.dose_per_fraction = np.hstack((self.dose_per_fraction, np.transpose([dose_1st])))

        # Fraction size update
        self.fractions.append(fractions)

        # Shuffle and repeat
        if fractions > 1:
            # Copy the 1st dose (N -1) times because the 1st dose is already
            # appended above
            d_group = np.vstack([dose_1st] * (fractions-1)).transpose()
            for i in range(0, fractions-1):
                # Shuffle to make fraction dose
                np.random.shuffle(d_group[0:,i])
            # Update the `dose_per_fraction` table
            self.dose_per_fraction = np.hstack((self.dose_per_fraction, d_group))

        # Update the final dose
        self.dose += np.sum(self.dose_per_fraction[0:, -fractions:], axis=1)

    def __add_dose(self, beam, compartment_id, start_time=0.0):
        """
        Beam may have multiple DVHs so add to each.

        Parameters
        ----------
        beam : tDVH
            The time-dependent DVH for the applied beams.
        compartment_ids : int
            The compartment ids to add the dose.
        start_time : float/int, optional
            The start time for the applied dose.

        Returns
        -------
        d : numpy.ndarray
            The applied doses to the blood particles.

        """
        d = np.zeros(self.blood_path.shape[0])
        assert(start_time + beam.dvh[0][-1] <= self.blood_path.shape[1] * self.dT)

        # For the index, time, and DVH function
        for i,(_,dvh) in enumerate(zip(beam.dvh[0], beam.dvh[1])):
            if dvh == None:
                continue
            # Two time-points for blood path (t0,t1) get irradiated by dose.
            t0 = start_time if i == 0 else start_time + beam.dvh[0][i-1]
            t1 = start_time + beam.dvh[0][i]
            dt1_t0 = t1 - t0
            # Calculate index of blood-particle distribution (BPD)
            lower_bound = np.floor(t0/self.dT)
            upper_bound = np.ceil(t1/self.dT)
            idx0 = int(lower_bound)
            idx1 = int(upper_bound)
            distance = idx1 - idx0

            for j, r in enumerate(self.blood_path.iloc[0:, idx0:idx1].values):
                # The 1st column: add dose before loop to account partial
                # delivery
                if r[0] == compartment_id:
                    d[j] += (self.dT - (t0 - lower_bound)) * dvh(np.random.uniform()) / dt1_t0
                # The path from second to last
                for _ in range(np.count_nonzero(r[1:] == compartment_id)):
                    d[j] += self.dT*dvh(np.random.uniform())/dt1_t0
                # The last column: minus dose in case of partial delivery
                if r[distance-1] == compartment_id:
                    d[j] -= (upper_bound - t1) * dvh(np.random.uniform()) / dt1_t0

        return d
