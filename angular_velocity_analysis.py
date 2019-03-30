
'''
This code has been written to examine what is the variation in angular velocity of orbit frame wrt inertial frame in orbit frame
'''

import numpy as np
import matplotlib.pyplot as plt


#declaring constants needed in this file
G = 6.67408e-11; #universal gravitational constant, SI
M_EARTH = 5.972e24; #mass of earth, kg

sgp_data = np.genfromtxt('sgp_output.csv',delimiter=",") #extracting positon data corresponding to the orbit for which we need to examine the variation
v_pos = sgp_data[:,1:4] #extracting the position values from sgp_data

v_pos_magnitude=np.zeros(120000)
v_w_OI_o=np.zeros((120000,2))

for i in range(0,120000):
	v_pos_magnitude[i] = np.linalg.norm(sgp_data[i,1:4]) #getting the value of distance of satellite from center earth
	v_w_OI_o[i,0] = np.sqrt(G*M_EARTH/(v_pos_magnitude[i])**3) #angular velocity of orbit frame wrt inertial frame in orbit frame

ang_vel_avg = np.average(v_w_OI_o[:,0])

for i in range(0,120000):
	v_w_OI_o[i,1] = (v_w_OI_o[i,0]-ang_vel_avg)/ang_vel_avg


plt.plot(sgp_data[0:120000,0],v_pos_magnitude)
plt.xlabel('time (in sec)' )
plt.ylabel('distance of satellite from center of earth (in m)')
plt.title('distance of satellite from center of earth versus time')
plt.show()

plt.plot(sgp_data[0:120000,0],v_w_OI_o[:,0])
plt.xlabel('time (in sec)' )
plt.ylabel('angular velocity of orbit frame wrt inertial frame in orbit frame (v_w_OI_o) (in rad/s)')
plt.title('v_w_OI_o versus time')
plt.show()

plt.plot(sgp_data[0:120000,0],v_w_OI_o[:,1])
plt.xlabel('time (in sec)' )
plt.ylabel('fractional error in angular velocity (v_w_OI_o) (in rad/s)')
plt.title('fractional error in v_w_OI_o versus time')
plt.show()