import numpy as np
import math


class Weibull:
    """
    Weibull distribution to calculate flow-out probability based on residence time.
    Clearly, one can readily implement other distributions by changing the pdf, cdf, etc. accordingly.
    Perhaps using scipy (which has all these functions implemented),
    (though it's typically slower compared to explicit implementation using Numpy).

    Note that with kappa = 1 (the shape parameter) the Weibull function reduces to the Exponential function.
    """

    def __init__(self, name, mtt, shape=2):
        """
        name  : name of the compartment (not used but nice to know).
        mtt   : mean transit time of the compartment (organ).
        shape : shape parameter of weibull distribution (k_parameter (k=1 corresponds to exponential))
        """
        self.name = name
        self.mtt = mtt
        # NOTE: the shape factor does not change the MTTs
        # (this is determined by the scales, see below). It only changes the shape of the transit time distribution.
        self.shape = shape
        # NOTE: the scales determine the speed at which a BP transits a specific organ
        # (by relating it to MTT so that the mean of the Weibull pdf equals the estimated MTT).
        self.scale = mtt / math.gamma(1 + 1 / self.shape)

    def pdf(self, t):
        """
        t : 1-array like
        PDF of Weibull distribution
        return probability of residence time
        """
        return (self.shape / self.scale) * (t / self.scale)**(self.shape - 1) * np.exp(-(t / self.scale)**self.shape)

    def cdf(self, t):
        """
        t : 1-d array like
        CDF of Weibull distribution
        return cumulative probability of residence time
        """
        return 1.0 - np.exp(-(t / self.scale)**self.shape)

    def cdf_inv(self, p):
        """
        p : 1-d array like
        Inverse of cumulative Weibull distribution
        return residence time at cumulative probability at p
        """
        return self.scale * pow(-1.0 * np.log(1.0 - p), 1.0/self.shape)

    def sf(self, t):
        """
        survival function.
        """
        return 1 - self.cdf(t)

    def hf(self, t):
        """
        hazard function.
        """
        return (self.shape / self.scale) * (t / self.scale)**(self.shape - 1)

    def initial_time_distribution(self, indices, option=1):
        """
        Initialize current dwell times in compartment based on the survival function.
        Two equivalent ways:
        1) drawing from normalized survival function
        2) draw of fractions of transit times weighted by the length of the interval
        (more chance for a particle to be in a long vs short interval).
        """

        if option == 1:
            max_val = self.cdf_inv(0.99)
            values = np.linspace(0, max_val, 1000)
            survival_probs = self.sf(values) / np.sum(self.sf(values))
            initial_times = np.random.choice(values, indices.size, replace=True, p=survival_probs)
        else:
            possible_tt = self.cdf_inv(np.linspace(1e-9, 1-1e-9, 1000))
            available_time = possible_tt / np.sum(possible_tt)
            transit_times = np.random.choice(possible_tt, indices.size, replace=True, p=available_time)
            initial_times = np.random.uniform(0, 1, indices.size) * transit_times
        return initial_times

    def is_leaving(self, t, dt):
        """
        t : an array of residence times.
        returns the indices of the elements that moving on to the next component

        What is the probability of transiting during this time step?
        It is the conditional probability that a BP hasn't transited thus far
        but will do so in the current time step; i.e. P(t <= T < t+dt | T>=t).
        This is given by the hazard function h(t) = pdf(t)/sf(t) multiplied by the time step.
        """
        random_numbers = np.random.uniform(size=t.size)
        # The following are equivalent (for small dt), but single call to hf faster:
        # p = (self.cdf(t + dt) - self.cdf(t)) / self.sf(t)
        p = self.hf(t) * dt
        return np.where(p > random_numbers)[0]

