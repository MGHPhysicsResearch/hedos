import copy

import numpy as np

from blooddvh import Weibull


class TimeDependentMarkovChain:
    __slots__ = ['prob', 'comp', 'size']
    def __init__(self, prob, scales, shapes=[]):
        """
        prob : numpy.ndarray
            The transition probabilities between compartments (organs).
        scales : numpy.ndarray
            The scaling factor for the Weibull distribution for each
            compartment (organ).
        shapes : numpy.ndarray, optional
            The shape factor for the Weibull distribution for each
            compartment (organ).

        Returns
        -------
        N/A

        """
        self.prob = copy.deepcopy(prob)
        self.size = prob.shape[0]
        # If `k`` is not given, `k`` is set to 1 for all compartments
        if len(shapes) == 0:
            shapes = np.ones(self.size)

        assert(len(scales) == self.size)

        # Array for Weibull distribution to determine leave or stay
        self.comp = []

        for row in range(self.size):
            # Create leaving (staying) probability distribution for a compartment
            self.comp.append(Weibull(scales[row], shapes[row], 1.0 - self.prob[row,row]))
            # Re-normalize the matrix
            self.prob[row, row] = 0.0
            self.prob[row] /= sum(self.prob[row])

    def walk(self, n_steps, c, t=0.0, dt=1.0):
        """
        Random walk for `n_steps` starting at c-th compartment.

        Parameters
        ----------
        c : int
            The compartment id to begin the random walk.
        t : float/int
            The time of initial aging at the compartment.
        dt : float/int
            The time step for the random walk.

        Returns
        -------
        steps : numpy.ndarray
            The simulated compartment ids for the random walk.

        """
        # Initialize the steps
        steps = np.zeros(n_steps, dtype=np.uint8)

        for i in range(n_steps):
            # Proceed step-by-step
            if self.comp[c].is_leaving(t):
                # Move to next compartment
                c = np.random.choice(self.size, p = self.prob[c])
                t = 0
            else:
                t += dt
            steps[i] = np.uint8(c)

        return steps
