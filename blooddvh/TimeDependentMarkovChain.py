import numpy as np
import copy 
from blooddvh import Weibull

class TimeDependentMarkovChain:
    __slots__ = ["prob", "comp", "size"]
    def __init__(self, prob, scales, shapes=[]):
        """
        prob   : 2d square matrix 
        scales : 1d array
        shapes : 1d array
        """
        self.prob = copy.deepcopy(prob)

        self.size = prob.shape[0]

        # if k is not given, k is set to 1 for all compartment
        if len(shapes) == 0 :
            shapes = np.ones(self.size)

        assert( len(scales) == self.size )

        # array for weibull distribution to determine leave or stay
        self.comp = []
        
        for row in range(self.size):
            
            # create leaving (staying) probability distribution for a compartment 
            self.comp.append( Weibull(scales[row], shapes[row], 1.0 - self.prob[row,row]) )

            # re-normalize matrix for 
            self.prob[row, row] = 0.0
            self.prob[row] /= sum(self.prob[row])

    def walk(self, n_steps, c, t=0.0, dt=1.0):
        """
        random walk for n_steps starting at c-th compartment
        the time t is initial aging at the compartment
        """

        # initialize steps
        steps = np.zeros(n_steps)
        
        for i in range(n_steps):
            # proceed step-by-step
            if self.comp[c].is_leaving(t):
                # move to next compartment
                c = np.random.choice(self.size, p = self.prob[c])
                t = 0
            else:
                t += dt
            
            steps[i] = c
        return steps
            
    
