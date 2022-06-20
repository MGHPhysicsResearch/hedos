import copy

import numpy as np


class Weibull:
    """
    Weibull distribution to calculate flow-out probability based on residence
    time.

    How to use:
    ```
    from blooddvh import Weibull

    liver = Weibull(18.0, k=2, p=0.05)  # Create a liver weibull model

    s = np.random.uniform(0,1, 10w000)  # Probability distribution
    x = liver.cdf_inv(s)                # Calculate age distribution in liver

    hist(x, bins=range(0,100), label=r'$mean \theta:{:.2f}$'.format(x.mean()))

    # Which blood particles are leaving the liver based on age
    b = liver.is_leaving(x)
    b = [int(liver.is_leaving(c)) for c in x]

    figure()
    hist(b)
    ```
    """
    __slots__ = ['scale', 'shape', 'prob']
    def __init__(self, l, k=1, p=0.5):
        """

        Parameters
        ----------
        l : float/int
            The scale parameter of the Weibull distribution (lambda).
        k : float/int, optional
            The shape parameter of the Weibull distribution (k = 1 for exponential).
        p : float/int, optional
            The probability of leaving.

        Returns
        -------
        N/A

        """
        self.scale = copy.deepcopy(l)
        self.shape = copy.deepcopy(k)
        self.prob = copy.deepcopy(p)

    def pdf(self, t):
        """
        Get the PDF of the Weibull distribution.

        Parameters
        ----------
        t : numpy.ndarray
            The times to determine the probability of leaving.

        Returns
        -------
        N/A : numpy.ndarray
            The return probabilities for the residence times.

        """
        return (self.shape / self.scale) * (t / self.scale)**(self.shape - 1) * np.exp(-(t /self.scale)**self.shape)

    def cdf(self, t):
        """
        Get the CDF of the Weibull distribution.

        Parameters
        ----------
        t : numpy.ndarray
            The times to determine the probability of leaving.

        Returns
        -------
        N/A : numpy.ndarray
            The cumulative return probabilities for the residence times.

        """
        return 1.0 - np.exp(-(t/self.scale)**self.shape)

    def cdf_inv(self, p):
        """
        Get the inverse of the Weibull distribution.

        Parameters
        ----------
        p : numpy.ndarray
            The probabilities to determine the time of leaving.

        Returns
        -------
        N/A : numpy.ndarray
            The residence time at cumulative probability for `p`.

        """
        return self.scale * pow( -1.0 * np.log(1.0 - p), 1.0/self.shape)

    def is_leaving(self, t):
        """
        Determine when the blood particle is leaving.

        Parameters
        ----------
        t : numpy.ndarray
            The times to determine the probability of leaving.

        Returns
        -------
        N/A : bool
            Return `True` if a particle with residence time `t` leaves. Return
            `False` if a particle with residence time `t` stays.

        """
        p = self.cdf(t) * self.prob/0.5
        if p > 1.0 :
            p = 1.0
        return np.random.choice([True, False], p=[p, 1.0-p])
