from workflows import BloodDoseFromFields, BloodDoseFromDVH

# =================== Simulation parameters ==================== #
sample_size = 100000    # nr simulation particles
nr_steps = 2400         # nr time steps
dt = 0.05               # in seconds
weibull_shape = 2
generate_new = True
random_walk = False
accumulate = False

simulation_parameters = {
    'sample_size': sample_size,
    'nr_steps': nr_steps,
    'dt': dt,
    'weibull_shape': weibull_shape,
    'generate_new': generate_new,
    'random_walk': random_walk,
    'accumulate': accumulate
}
# ============================================================== #

# ==================== Patient parameters ====================== #
gender = 'F'
tumor_site = 'lung'
tumor_volume_fraction = 0.05
relative_blood_density = 1.0
relative_perfusion = 1.0
organs = ['lung', 'heart', 'aorta', 'pulmonary_artery', 'superior_vena_cava', 'tumor']
# organs = ['lung', 'heart', 'tumor']

print('Gender: ', gender)
if gender.lower() in ['m', 'male']:
    sheet_name = 'male'
    TBV = 5.3  # 5.3 L total simulation volume.
    CO = 6.5 / 60  # 6.5 L/min cardiac output.
elif gender.lower() in ['f', 'female']:
    sheet_name = 'female_new'
    TBV = 3.9  # 3.9 L total simulation volume.
    CO = 5.9 / 60  # 5.9 L/min cardiac output.
else:
    raise ValueError('Cannot deduce gender.')

patient_parameters = {
    'gender': gender,
    'sheet_name': sheet_name,
    'tumor_site': tumor_site,
    'tumor_volume_fraction': tumor_volume_fraction,
    'relative_blood_density': relative_blood_density,
    'relative_perfusion': relative_perfusion,
    'organs': organs,
    'TBV': TBV,
    'CO': CO
}
# ============================================================== #

# ==================== Treatment parameters ==================== #
nr_fractions = 30
# total_beam_on_time = 80
# start_times = [10, 40, 70, 100]
# beam_on_times = [20, 20, 20, 20]
total_beam_on_time = 10
start_times = [10]
beam_on_times = [10]
assert(sum(beam_on_times) == total_beam_on_time), 'Beam-on-time of separate fields should equal total beam-on-time.'
assert(start_times[:-1] + beam_on_times[:-1] <= start_times[1:]), 'Cannot start new field before completing current.'

treatment_parameters = {
    'nr_fractions': nr_fractions,
    'total_beam_on_time': total_beam_on_time,
    'start_times': start_times,
    'beam_on_times': beam_on_times,
}
# ============================================================== #

# BloodDoseFromFields.blood_dose_distribution(simulation_parameters, patient_parameters, treatment_parameters)
BloodDoseFromDVH.blood_dose_distribution(simulation_parameters, patient_parameters, treatment_parameters)
