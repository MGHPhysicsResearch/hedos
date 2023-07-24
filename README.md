<img src="figures/logo.png">
HEDOS : Hematological Dose

## Installation

Checkout the python source code
```bash
$ git clone https://github.com/mghro/hedos.git

```
Install dependent packages

```bash
$ pip3 install -r requirements.txt
```

## Background

The modules developed for this project were the following:
1. `CompartmentModel` reads in organ's transition flow and generate transition network.
2. `BloodDistribution` receives the Markov model from `CompartmentModel` and produces blood particle paths. Further, it can save the generated paths and read it back.
3. `tDVH` stores the DVH as a function of time.
4. `bDVH` computes blood dose by applying `tDVH` to `BloodDistribution`.

Their relationship looks like following image:
<img src="figures/blood_diagram_hedos.png">

## Getting started

These commands can be adjusted based off of the `tests/test.ipynb` file.

### 1. Load python module

```python
import os

import matplotlib.pyplot as plt

from blooddvh import CompartmentModel
from blooddvh import BloodDistribution
from blooddvh import tDVH
from blooddvh import bDVH
```

### 2. Create a compartment model from an the ICRP 89 Excel sheet
```python
sample_size = 100
time_per_step = 1   # Seconds
steps_per_min = 60  # Number of steps per min, e.g.,  step resolutions are 1 sec and 0.1 sec for 60 and 600, respectively
model = CompartmentModel(os.path.join('input', 'ICRP89_compartment_model.xlsx'), 'male', vol=5.3, cardiac=6.5, resolution=steps_per_min)
```

### 3. Generate blood distribution using Markov chain
```python
blood = BloodDistribution()
blood.generate_from_markov(model.markov, model.name, model.volume, time_per_step, sample_size, steps_per_min)
```

### 4. Create a time-dependent DVH

```python
dose = tDVH()
# First 10 sec, 2 Gy uniform
dose.add(10, lambda x: 2)
# Next 10 sec, no dose
dose.add(10, None)
# Next 10 sec, 5 Gy uniform
dose.add(10, lambda x: 5)
```

### 5. Apply tDVH to blood path
tDVH can be applied to multiple organ and different time points
```python
# To choose an organ where the dose is delivered, 0: brain, 19: liver
print(list(zip(range(len(model.name)), model.name)))
blood_dose = bDVH(blood.df, blood.dt)
# Add dose to Liver
blood_dose.add_dose(dose, 19, start_time=0)
blood_dose.add_dose(dose, 19, start_time=4)
# Add dose to Brain
blood_dose.add_dose(dose, 0, start_time=0)
# Draw your blood DVH (differential)
plt.figure(figsize=(10,6))
plt.hist(blood_dose.dose)
plt.xlabel('Dose (Gy)', fontsize=16)
plt.ylabel('# of Events', fontsize=16)
```

[(0, 'brain'), (1, 'thyroid'), (2, 'breast'), (3, 'lymph_nodes'), (4, 'large_veins'), (5, 'all_other'), (6, 'fat'), (7, 'skeletal_muscle'), (8, 'adrenals'), (9, 'skin'), (10, 'red_marrow'), (11, 'spongy_bone'), (12, 'compact_bone'), (13, 'other_skeleton'), (14, 'bronchial'), (15, 'pulmonary'), (16, 'right_heart'), (17, 'left_heart'), (18, 'coronary'), (19, 'liver'), (20, 'kidneys'), (21, 'bladder'), (22, 'gonads'), (23, 'aorta_large_arteries'), (24, 'pancreas'), (25, 'spleen'), (26, 'stomach_oesophagus'), (27, 'small_intestine'), (28, 'large_intestine')]

<img src="figures/getting_started_results.png">

### 6. Further information (see the [test notebook](tests/test.ipynb) for more)

## Team
- Chris Beekman
- Jungwook Shin
- Stella Xing
- Lucas McCullum
- Clemens Grassberger
- Harald Paganetti

## Acknowledgements
This work was supported by:
- R21 CA248118 : A Computational Method to Calculate the Radiation Dose to Circulating Lymphocytes
- R01 CA248901 : Developing whole-body computational phantoms for blood dosimetry to model the impact of radiation on the immune system

## Publication
Shin J, Xing S, McCullum L, et al. HEDOS-a computational tool to assess radiation dose to circulating blood cells during external beam radiotherapy based on whole-body blood flow simulations. Phys Med Biol. 2021;66(16):10.1088/1361-6560/ac16ea. Published 2021 Aug 3. doi:10.1088/1361-6560/ac16ea
