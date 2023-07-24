import numpy as np
from scipy import interpolate
import pandas as pd


def field_to_func(vol, gridpoints):
    """
    Create function from volume values
    Parameters
    ----------
    vol: ndarray (3D), required.
        volume to be used for the interpolation.
    gridpoints: tuple of x, y, z coordinates, required.

    Returns
    -------
    field_fn: function.
        function representation of the volume.
    """
    field_fn = interpolate.RegularGridInterpolator(gridpoints, vol, bounds_error=False, fill_value=0)
    return field_fn


class DoseRate:
    """
    Get the dose rate, either as a dose rate histogram, or as an interpolated function over the segmentation volume.
    """
    def __init__(self, patient, n_fractions, total_beam_on_time):
        self.patient = patient
        self.n_fractions = n_fractions
        self.total_beam_on_time = total_beam_on_time
        self.dose_rate = None
        self.dose_rate_func = None
        self.MBD = None

    def get_mean_organ_dose(self, organ_name, accumulate=False):
        if accumulate:
            return self.patient.get_mean_organ_dose(organ_name)
        else:
            return self.patient.get_mean_organ_dose(organ_name) / self.n_fractions

    def get_dose_rate(self):
        self.dose_rate = self.patient.dose / self.n_fractions / self.total_beam_on_time
        self.dose_rate_func = field_to_func(self.dose_rate, self.patient.gridpoints)

    def get_dose_rate_hist(self, organ_name):
        bin_size = 0.1 / self.n_fractions / self.total_beam_on_time
        bin_max = np.ceil(np.max(self.patient.dose)) / self.n_fractions / self.total_beam_on_time
        bins = np.arange(0, bin_max + bin_size, bin_size)
        idx = np.where(self.patient.seg_organs[organ_name] == 1)
        dose_rate = self.dose_rate[idx[0], idx[1], idx[2]]
        frequency, dose_rate_bins = np.histogram(dose_rate, bins=bins, density=True)
        frequency *= bin_size * 100
        return frequency, dose_rate_bins

    def calculate_mean_blood_dose(self, total_blood_volume, blood_volumes, accumulate=False):
        mean_organ_doses = np.array([self.get_mean_organ_dose(organ_name, accumulate=accumulate)
                                     for organ_name in self.patient.seg_organs.keys()])
        self.MBD = np.sum(blood_volumes * mean_organ_doses) / total_blood_volume


class DoseRateFromDVH:
    """
    Instead, use DVH files directly. This could be advantageous for data-sharing between institutes.
    """
    def __init__(self, n_fractions, total_beam_on_time):
        self.n_fractions = n_fractions
        self.total_beam_on_time = total_beam_on_time

    def get_dose_rate_hist(self, dvh_file):
        dvh = pd.read_csv(dvh_file, header=0)
        dose_rate_bins = dvh['dose_bins'].values / self.n_fractions / self.total_beam_on_time
        coverage = dvh['coverage'].values
        frequency = np.diff(coverage[::-1])[::-1]
        return frequency, dose_rate_bins
