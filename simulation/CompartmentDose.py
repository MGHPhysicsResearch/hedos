import copy
import numpy as np
import numbers
from scipy.sparse import csc_matrix
from scipy.spatial import KDTree


class CompartmentDose:
    def __init__(self, blood_path, dt):
        """
        blood_path: ndarray of the spatiotemporal distribution of simulation particles
        dt        : Time (sec) per step in blood_path
        """
        self.blood_path = blood_path
        self.dt = dt

        # total dose
        self.dose = np.zeros(self.blood_path.shape[0])

    def volume_gt_dose(self, threshold):
        """
        threshold: threshold dose (Gy)
        return volume (%) of BP g.t threshold, mean dose, std dose
        """
        # return mean/std of dose?
        dose_above_threshold = self.dose[self.dose > threshold]
        volume_above_threshold = 100.0 * dose_above_threshold.size / self.dose.size
        return volume_above_threshold, dose_above_threshold.mean(), dose_above_threshold.std()

    def dose_at_top_volume(self, top_n):
        """
        top_n (%) : upper volume, e.g., 2 for 2 %
        return dose (Gy) at volume, mean dose, std dose
        """
        top_percentile = np.percentile(self.dose, q=100 - top_n)
        top_n_dose = self.dose[self.dose > top_percentile]
        return top_percentile, top_n_dose.mean(), top_n_dose.std()

    def add_dose(self, dose_rate, compartment_ids, start_time, beam_on_time):
        """
        dose_rate           : either a value (homogeneous dose) or dose histogram (heterogeneous dose distribution).
        compartment_ids     : list of compartment ids that gets the given dose; the idea that this could be more than
                              one is for instance when you have multiple substructures,
                              but you only have the DVH for their union.
        start_time          : time-point beam starts to delivery
        beam_on_time        : duration of the applied dose
        """
        assert (isinstance(dose_rate, numbers.Number) or isinstance(dose_rate, tuple)), \
            'dose_function needs to be either a value (homogeneous dose) ' \
            'or the output of np.histogram (a tuple for heterogeneous dose).'
        assert (start_time + beam_on_time <= self.blood_path.shape[1] * self.dt)
        if compartment_ids is None:
            return 0
        if not isinstance(compartment_ids, list):
            compartment_ids = [compartment_ids]

        # Calculate index of simulation-particle distribution (BPD)
        idx0 = int(np.floor(start_time / self.dt))
        idx1 = int(np.ceil((start_time + beam_on_time) / self.dt))

        # don't copy the entire simulation path; make use of the sparsity to reduce memory burden:
        # This gives sparse representation of (n_particles x n_timesteps),
        # indicating for each particle whether in compartment or not at that time.
        in_compartment = csc_matrix(sum([self.blood_path[:, idx0:idx1] == compartment_id
                                         for compartment_id in compartment_ids]))

        if isinstance(dose_rate, numbers.Number):
            d_dose = self.dt * dose_rate * np.ones(shape=in_compartment.data.size)
        else:
            dose_values = (dose_rate[1][1:] + dose_rate[1][:-1] - np.diff(dose_rate[1])) / 2
            d_dose = self.dt * np.random.choice(dose_values, size=in_compartment.data.size,
                                                replace=True, p=dose_rate[0]/np.sum(dose_rate[0]))
        # Idea: you work with a "flattened array", but the csc-matrix remembers which particle is which...
        in_compartment.data = in_compartment.data * d_dose
        # ... so that we can sum over it (dose-accumulation for each particle).
        self.dose += np.array(in_compartment.sum(axis=1)).flatten()

    def add_dose_random_walk(self, dose_rate_func, compartment_ids, start_time, beam_on_time):
        """
        This function accumulates dose for simulation particles appearing and disappearing into and out of the compartment.
        Once they appear, they embark on a random walk through region confined by a supplied segmentation.
        Make sure to call the 'prepare' function first to define the segmentation and the KDTree that is used here.

        dose_rate_func : function of dose rate.
        step_size : typical step size of the motion of the BP performing a random walk
        compartment_ids : list of compartment ids that gets the given dose
        start_time : time-point beam starts to delivery
        beam_on_time : duration of the applied dose
        """
        assert (start_time + beam_on_time <= self.blood_path.shape[1] * self.dt)
        if compartment_ids is None:
            return 0
        if not isinstance(compartment_ids, list):
            compartment_ids = [compartment_ids]

        # Calculate index of simulation-particle distribution (BPD)
        idx0 = int(np.floor(start_time / self.dt))
        idx1 = int(np.ceil((start_time + beam_on_time) / self.dt))
        nr_time_steps = idx1 - idx0

        # initialize:
        pos = np.full(shape=(self.blood_path.shape[0], 3), fill_value=np.nan)
        all_idx = np.arange(self.positions.shape[0])
        indices_old = np.array([], dtype=int)

        in_compartment = np.sum([self.blood_path[:, idx0:idx1] == compartment_id for compartment_id in compartment_ids],
                                axis=0, dtype=bool)

        accept_perc = []
        for step in range(nr_time_steps):
            indices = np.where(in_compartment[:, step])[0]
            indices_inject = np.setdiff1d(indices, indices_old, assume_unique=True)
            indices_eject = np.setdiff1d(indices_old, indices, assume_unique=True)
            idx = np.random.choice(all_idx, size=indices_inject.size)
            pos[indices_inject] = self.positions[idx]
            # this is of course not entirely correct, you would want to sample the step_size,
            # uniformly picking the direction. Good enough though (it's an approx anyway):
            dist = np.random.uniform(0, self.step_size, size=indices_old.size * 3).reshape(indices_old.size, 3)
            # check if particles are still inside organ, reject the moves for those who landed outside.
            d = self.kd_tree.query(pos[indices_old] + dist, k=1)[0]
            idx_accept = np.where(d < self.d_max)[0]
            if d.size > 0:
                accept_perc.append(idx_accept.size / d.size)
            # move the particles:
            pos[indices_old[idx_accept]] += dist[idx_accept]
            # discard the particles that have left the compartment:
            pos[indices_eject] = np.nan
            indices_old = indices
            # accumulate dose
            self.dose[indices] += self.dt * dose_rate_func(pos[indices])
        print('Percentage accepted: {:.2f}%'.format(sum(accept_perc)/len(accept_perc) * 100))
        print('Dose added.')

    def _prepare(self, grid, seg, down_sample):
        """
        Prepare for the random walk.
        Velocity is a free parameter whose influence can be studied,
        higher values lead to bigger jumps with every time step,
        i.e. velocity=0 --> stationary simulation particles while in compartment.
        Of course, in reality, simulation flow has much more directionality than just random,
        could potentially be implemented with Levy flights or by adding momentum...
        """
        if down_sample is not None:
            grid = tuple(grid[i][::down_sample[i]] for i in range(3))
            seg = seg[::down_sample[0], ::down_sample[1], ::down_sample[2]]
        seg_idx = np.where(seg == 1)
        gridpts = np.stack(np.meshgrid(*grid, indexing='ij'), axis=-1)
        self.positions = gridpts[seg_idx[0], seg_idx[1], seg_idx[2]]
        self.kd_tree = KDTree(self.positions)
        # tolerance:
        gridspacing = np.array([np.diff(grid[i])[0] for i in range(3)])
        self.d_max = np.sqrt(np.sum(np.square(0.5 * gridspacing)))
        # velocity of particles:
        velocity = 20
        self.step_size = velocity * self.dt

    def prepare(self, grid, seg, down_sample=None):
        self._prepare(grid, seg, down_sample=down_sample)

    def repeat(self, n_fractions):
        """
        Dose accumulation over multiple fractions.
        With every fraction, the dose array is shuffled and added to itself (thus assuming total mixing).
        For n_fractions large, the resulting dose distribution will become normal (CLT).
        """
        f_dose = copy.deepcopy(self.dose)
        self.dose = []
        for _ in range(n_fractions):
            np.random.shuffle(f_dose)
            self.dose.append(copy.deepcopy(f_dose))
        self.dose = sum(self.dose)
