from scipy.io import loadmat

# In Abdel's code, three data were given.
# 'fixedData' is for CT but never used
# 'Dose_field' is for Dose
# 'Trajectory' is for trajectory points
# 'movingData' is for registration output (not used)

# This Dose_field.mat has three doses
# three doses
# dose[0][0] - dose[0][2]
# shape(dose[0][0]) (200, 512, 512)
# time_delivery = 60 sec
# fractionN = 1
# dose =
dose = loadmat("../input/Dose_field.mat")['Dose_field'][0]
shape(dose)[0]
delivery_time = 60 
dose /= delivery_time #dose-rate

###dose[0] 

ct   = loadmat("../input/fixedData.mat")['fixedData'][0][0][0]
#ct[0][0][0] hu map
# shape(dose_header[0][0][0])  (200, 512, 512)
# resolution? ct[ array([[0.684],  [0.684]]) check with Abdel


moving = loadmat("../input/movingData.mat")['movingData']
# 3D image with dicom tags
# shape(160, 256, 212)


# these are dictionaries
# (1050,2)
# row: trajectories 1050 trajectories,
# column: 0 -> dummy, 1-> points
# Q: connectivity? branch?
# shape(traj) : 1050 tracks
# Trajectory_V19CT.mat wasn't used.
traj = loadmat("../input/Trajectory.mat")['Trajectory']
# id : 0-1050
# 
# traj[0][1][0:10,0] : first trajetory, 1 is dummy, 10 x points

#
nb_traj = shape(traj)[0]
nb_pts_traj = np.array([ len(traj[p][1]) for p in range(nb_traj) ] )


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = figure()
ax = fig.add_subplot(111, projection='3d')
t_id = 100 ;
for i in range( 0, nb_traj) : 
    ax.plot(traj[i][1][0:,0], traj[i][1][0:,1], traj[i][1][0:,2])


# Blood vol and flow 
# matlab function : Blood_vol2flowV6

#TrackID_spacingV6.m

# repmat(resCT , [size(path,1),1]); #repeat copy of matrix
