import numpy as np

from simulation import ExpandFlowModel, TemporalDistribution, DoseRateFromDVH, CompartmentDose
from PlotDoseDistribution import plot_dose_distribution


def blood_dose_distribution(simulation_params, patient_params, treatment_params):

    # ======= Step 1. Initialize simulation ============================= #
    filename = '../input/phantom/ICRP89_compartment_model.xlsx'
    model = ExpandFlowModel(filename, patient_params, simulation_params)
    # add the tumor:
    if 'tumor' in patient_params['organs']:
        model.split_box_parallel('tumor', box_dict=patient_params)
    # ============================================================== #

    # ======== Step 2. Generate distribution ======================= #
    blood = TemporalDistribution(model)
    if simulation_params['generate_new']:
        model.construct_weibull()
        blood.generate_from_weibull()
        blood.save('../input/blood_path.npy')

        # Could also do a Markov process, i.e. corresponding to exponential transit time distribution
        # This is the same as the above with Weibull shape_parameter=1.
        # model.construct_markov()
        # blood.generate_from_markov()
    else:
        blood.load('../input/blood_path.npy')

    # blood.plot_time_distributions(['lung'])
    # blood.plot_inflow_outflow(['lung'])
    # blood.plot_volumes_over_time()
    # blood.plot_normalized_volume_over_time(['lung', 'tumor'])
    # blood.plot_final_blood_volumes()
    # ============================================================== #

    # ======== Step 3. Accumulate dose ============================= #
    dose = DoseRateFromDVH(n_fractions=treatment_params['nr_fractions'],
                           total_beam_on_time=treatment_params['total_beam_on_time'])

    compartment_ids = [[i for i, name in enumerate(model.names) if organ in name] for organ in patient_params['organs']]
    dose_contributions = {}
    for organ, compartment_id in zip(patient_params['organs'], compartment_ids):
        blood_dose = CompartmentDose(blood.path, model.dt)
        dose_rate_hist = dose.get_dose_rate_hist('../input/patient/DVHs/' + organ + '_DVH.csv')
        for start_time, beam_on_time in zip(treatment_params['start_times'], treatment_params['beam_on_times']):
            blood_dose.add_dose(dose_rate_hist, compartment_id, start_time=start_time, beam_on_time=beam_on_time)
        dose_contributions[organ] = blood_dose.dose

    blood_dose_total = CompartmentDose(blood.path, model.dt)
    blood_dose_total.dose = sum(list(dose_contributions.values()))
    if simulation_params['accumulate']:
        blood_dose_total.repeat(treatment_params['nr_fractions'])

    plot_dose_distribution(blood_dose_total, dose_contributions)
    # ============================================================== #

    from scipy.stats import wasserstein_distance

    max_val = np.max(blood_dose_total.dose)
    step_size = max_val / 1000
    weights, bins = np.histogram(blood_dose_total.dose,
                                 bins=np.linspace(0, max_val+step_size, 1002),
                                 density=True)
    values = bins[:-1]
    weights_ref = np.zeros_like(weights)
    weights_ref[0] = 1

    w = wasserstein_distance(values, values, weights_ref, weights)
    w2 = np.sum((1 - np.cumsum(weights * step_size)) * step_size)
    print('Wasserstein distance: {:.3f}'.format(w))
    print('Wasserstein distance: {:.3f}'.format(w2))
    print('Mean dose: {:.3f}'.format(np.mean(blood_dose_total.dose)))
