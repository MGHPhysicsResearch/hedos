from blooddvh import CompartmentModel, save_blood_distribution

sample_size    =   10000 #10000
steps_per_min  =      60 # per min, e.g.,  step resolutions are 1 sec and 0.1 sec for 60 and 600, respectively

# 1. Load your flow transition map from an excel file and tab name
# total blood volume : 5.3 L
# cardiac output: 6.5 L/min
# time-resolution is 60 per min, i.e., 1 sec
model = CompartmentModel("../input/ICRP89_compartment_model.xlsx", "male", vol=5.3, cardiac=6.5, resolution=steps_per_min)

# 2. Create BloodDistribution
blood = BloodDistribution() #

# option 1: create from markov chain
blood.generate_from_markov(model.markov, model.name, model.volume, 1.0, 1000, steps_per_min)

# option 2: load from pre-generated 
#blood.read_from_excel("BP_10K_60_30mins.xlsx", "BP")
# 3.2 save blood distribution to file
#blood.save_blood_distribution(blood_df, "BP_10K_60_30min.xlsx")


# 3. Load DVH(s) with time-information
