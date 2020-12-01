import numpy as np
import copy
import math 

class Weibull:
    """
    Weibull distribution to calculate flow-out probability based on residence time.

    How to use:
    
    from blooddvh import Weibull

    liver = Weibull(18.0, k=2, p=0.05)   # create a liver weibull model
    
    s = np.random.uniform(0,1, 10w000)    # probability distribution
    x = liver.cdf_inv(s)                 # calculate age distribution in liver

    hist(x, bins=range(0,100), label=r'$mean \theta:{:.2f}$'.format(x.mean()))
    
    b = liver.is_leaving(x)  # who is leaving liver based on age
    b = [ int(liver.is_leaving(c)) for c in x]  # 

    figure()
    hist(b)
    """
    __slots__ = ["scale", "shape", "prob"]
    def __init__(self, l, k=1, p=0.5):
        """
        scale : scale parameter of weibull distribution (lambda) 
        shape : shape parameter of weibull distribution (k  and 1 for exponential)
        prob  : leaving probability by default 0.5
        """
        self.scale = copy.deepcopy(l)
        self.shape = copy.deepcopy(k)
        self.prob  = copy.deepcopy(p)

    def pdf(self, t):
        """
        t : 1-array like
        PDF of Weibull distribution 
        return probability of residence time
        """
        return (self.shape / self.scale) * (t / self.scale)**(self.shape - 1) * np.exp(-(t /self.scale)**self.shape)

    def cdf(self, t):
        """
        t : 1-d array like
        CDF of Weibull distribution
        return cumulative probability of residence time
        """
        return 1.0 - np.exp(-(t/self.scale)**self.shape)

    def cdf_inv(self, p):
        """
        p : 1-d array like
        Inverse of Weibull distribution
        return residence time at cumulative probability at p
        """
        return self.scale * pow( -1.0 * np.log(1.0 - p), 1.0/self.shape)


    def is_leaving(self, t):
        """
        t : "a value" not array
        return true if a particle with residence time t leaves
        return false if a particle with resience time t stays
        """
        p = self.cdf(t) * self.prob/0.5
        if p > 1.0 :
            p = 1.0
        return np.random.choice([True, False], p=[p, 1.0-p]) 
