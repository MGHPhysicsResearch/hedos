import numpy as np
import copy

from simulation import Weibull


class Chain:
    def __init__(self, names, prob, mtt, dt, k):
        """
        names  : list of names of compartments
        prob   : 2d square matrix
        scales : 1d array
        shapes : 1d array
        """
        # some checks:
        assert(prob.shape[0] == prob.shape[1] == len(names) == mtt.size), 'Dimensions do not match.'
        self.prob = copy.deepcopy(prob)
        self.size = prob.shape[0]
        self.dt = dt
        self.progress = 0

        # list of weibull distributions to determine leave or stay
        self.comp = [Weibull(names[row], mtt[row], shape=k) for row in range(self.size)]

        # re-normalize matrix to determine where to go, given a particle leaves the compartment
        np.fill_diagonal(self.prob, 0.0)
        self.prob /= np.sum(self.prob, axis=1, keepdims=True)

        # list of tuples(compartment, probability) of successors for each compartment:
        self.options = [(np.arange(self.size)[p > 0], p[p > 0]) for p in self.prob]

    def walk_v1(self, n_steps, c):
        """
        Walk through time and compartments of all particles in simulation.
        The time t is initial aging at the compartment.
        It is initialized so that system is in equilibrium right from the start.
        """
        path = np.empty(shape=(c.size, n_steps+1), dtype=np.uint8)
        path[:, 0] = c

        # initialize dwell times of particles:
        t = np.empty(shape=c.size, dtype=np.float32)
        for i, comp in enumerate(self.comp):
            indices = np.where(c == i)[0]
            t[indices] = comp.initial_time_distribution(indices)

        # loop over time steps:
        for step in range(n_steps):
            # for each compartment, determine which particles leave and where they go to:
            for i, comp in enumerate(self.comp):
                indices = np.where(c == i)[0]
                change_indices = indices[comp.is_leaving(t[indices], dt=self.dt)]
                t[indices] += self.dt
                if change_indices.size > 0:
                    if self.options[i][0].size > 1:
                        c[change_indices] = np.random.choice(self.options[i][0], size=change_indices.size,
                                                             p=self.options[i][1])
                    else:
                        c[change_indices] = self.options[i][0][0]
                    t[change_indices] = 0
            path[:, step + 1] = np.uint8(c)
            # print progress:
            self._print_progress(n_steps, step)
        return path

    def walk_v2(self, n_steps, c):
        """
        random walk for n_steps starting at c-th compartment
        the time t is initial aging at the compartment
        This produces the same thing as above but differently.
        It compares the current dwell time with a particles transit time,
        which is determined upon entrance of the compartment.

        NOTE: watch out with indexing; don't do something like t[indices1][indices2] = arr
        (the first indexing creates a copy instead of a view). Instead, do t[indices1[indices2]] = arr.
        """
        path = np.empty(shape=(c.size, n_steps+1), dtype=np.uint8)
        path[:, 0] = c

        # initialize transit times and initial times:
        tt = np.empty(shape=c.size, dtype=np.float32)
        t = np.empty(shape=c.size, dtype=np.float32)
        for i, comp in enumerate(self.comp):
            indices = np.where(c == i)[0]
            possible_tt = comp.cdf_inv(np.linspace(1e-9, 1-1e-9, 1000))
            available_time = possible_tt / np.sum(possible_tt)
            tt[indices] = np.random.choice(possible_tt, indices.size, replace=True, p=available_time)
            t[indices] = np.random.uniform(0, 1, indices.size) * tt[indices]

        for step in range(n_steps):
            for i, comp in enumerate(self.comp):
                indices = np.where(c == i)[0]
                new_indices = indices[np.where(t[indices] == 0)[0]]
                tt[new_indices] = comp.cdf_inv(np.random.uniform(size=new_indices.size))
                change_indices = indices[np.where(t[indices] > tt[indices])[0]]
                t[indices] += self.dt
                if change_indices.size > 0:
                    if self.options[i][0].size > 1:
                        c[change_indices] = np.random.choice(self.options[i][0], size=change_indices.size,
                                                             p=self.options[i][1])
                    else:
                        c[change_indices] = self.options[i][0][0]
                    t[change_indices] = 0
            path[:, step + 1] = np.uint8(c)
            # print progress:
            self._print_progress(n_steps, step)
        return path

    def _print_progress(self, n_steps, step):
        percentage_done = int(100 * step / (n_steps - 1))
        if not self.progress == percentage_done:
            if percentage_done % 10 == 0:
                print('{:}% of blood flow simulation done.'.format(int(percentage_done)))
                self.progress = percentage_done


class MarkovChain:
    def __init__(self, names, prob):
        """
        names  : list of names of compartments
        prob   : 2d square matrix
        """
        # some checks:
        assert (prob.shape[0] == prob.shape[1] == len(names)), 'Dimensions do not match.'
        self.names = names
        self.prob = copy.deepcopy(prob)
        self.size = prob.shape[0]

        # leaving probabilities:
        self.p_leaving = 1 - np.diag(self.prob)

        # re-normalize matrix to determine where to go, given a particle leaves the compartment
        np.fill_diagonal(self.prob, 0.0)
        self.prob /= np.sum(self.prob, axis=1, keepdims=True)

    def is_leaving(self, size, c_id, random_numbers):
        """
        size: how many particles?
        c_id: id of the compartment
        """
        p = np.ones(shape=size) * self.p_leaving[c_id]
        return p > random_numbers

    def walk(self, n_steps, c):
        path = np.zeros(shape=(c.size, n_steps+1), dtype=np.uint8)
        path[:, 0] = c

        for step in range(n_steps):
            c_copy = copy.deepcopy(c)
            random_numbers = np.random.uniform(0, 1, c.size)
            for i in range(self.size):
                indices = np.where(c == i)[0]
                change = self.is_leaving(indices.size, i, random_numbers[indices])
                c_copy[indices[change]] = np.random.choice(self.size, size=np.sum(change), p=self.prob[i])
            c = c_copy
            path[:, step + 1] = np.uint8(c)
        return path
