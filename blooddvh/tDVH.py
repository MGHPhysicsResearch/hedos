import numpy as np
from scipy import interpolate


"""
Time-dependent DVH class consisting of dT and DVH.
"""
class tDVH:
    __slots__ = ['dvh', 'seed', 'start_time']

    def __init__(self):
        """
        dvh[0] : time, e.g., elapsed time (array) from start
        dvh[1] : dvh function (list), if no dose for time, None is recommended instead of 0 due to performance issue

        For example,
        dvh = [ np.array([1,3,6,7,10]), [lambda v:3, lambda v:0, v2d, v2d_const, None]]

        """
        self.dvh = [np.array([]),[]]
        self.start_time = 0.0

    def add(self, beam_on_time, dvh_function):
        """
        Fill dvh functions as a function of time

        Parameters
        ----------
        beam_on_time : float/int
            The delivery time of this DVH in seconds. The time is relative to
            beam start.
        dvh_function : lambda or intpl1d
            How to apply the DVH as a function.

        Returns
        -------
        N/A

        """
        stop_time = 0.0 if self.dvh[0].size == 0 else self.dvh[0][-1]
        self.dvh[0] = np.append(self.dvh[0], stop_time+beam_on_time)
        self.dvh[1].append(dvh_function)

    def add_array(self, beam_on_time, dvh_points=np.array([]), max_value=100.0):
        """
        Parameters
        ----------
        beam_on_time : float/int
            The delivery time for the DVH.
        dvh_points : numpy.ndarray, optional
            2D-dose and probability where the x-axis is dose (Gy) and the
            y-axis is volume (0-100%).
        max_value : float/int
            The value to normalize volume of DVH to 0-1 scale.

        Returns
        -------
        N/A

        """
        if dvh_points.shape != (0,):
            v2d = interpolate.interp1d(dvh_points[0:,1]/max_value, dvh_points[0:,0], fill_value=(0.0, 0.0), bounds_error=False)
            self.add(beam_on_time, v2d)
        else:
            self.add(beam_on_time, None)
