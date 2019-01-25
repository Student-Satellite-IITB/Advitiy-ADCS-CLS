import numpy as np
import frames as fs
from constants_1U import R_EARTH
import math
import matplotlib.pyplot as plt

'''
	This code generates csv file for latitude, longitude, altitude (LLA) corresponding to time data.
	ECI frame position is provided from sgp_output.csv
		Time: (since epoch) in seconds
		latitude: -90 to 90 degrees
		longitude: -180 to 180 degrees (-180 excluded)
		altitude: in meters
'''

m_sgp_output_i = np.genfromtxt('sgp_output.csv', delimiter=",")
N = m_sgp_output_i.shape[0]
m_sgp_ecef = np.zeros([N,4])
m_LLA = np.zeros([N,4])

for k in range(0,N):
	if math.fmod(k,N/100) == 0:
		print (int(100.*k/N)) 

	v_i = m_sgp_output_i[k,1:4]		#iniertial frame position
	time = m_sgp_output_i[k,0]		#time in sec
	#get position in ecef 	
	v_ecef = fs.ecif2ecef(v_i,time)
	#get latitude and longitude and altitude
	v_latlon = fs.latlon(v_ecef.copy())
	alt = np.linalg.norm(v_i.copy()) - R_EARTH	#in meters
	
	m_sgp_ecef[k,0] = time
	m_LLA[k,0] = time
	m_sgp_ecef[k,1:4] = v_ecef
	m_LLA[k,1:4] = np.append(v_latlon,alt)
np.savetxt('LLA.csv',m_LLA, delimiter=",")
print("LLA done")