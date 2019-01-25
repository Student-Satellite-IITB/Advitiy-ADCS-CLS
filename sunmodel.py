import numpy as np
from constants_1U import EPOCH, EQUINOX
import scipy.io as sio
'''
	This code gives sun-vector (unit vector from center of earth to center of sun) in ECI frame.
	Input file is sgp data with first column as time series in seconds
	The report for this code is uploaded on thread : Controls | Environment model - EV_sunmodel-report_01.pdf
'''


m_sgp_output = np.genfromtxt('sgp_output.csv', delimiter=",")

T = m_sgp_output[:,0].copy()	#time in seconds
N = len(T)
m_si_output = np.zeros((N,4))

#The time from equinox and the initial time (launch) in seconds
initialdelay = (EPOCH - EQUINOX).total_seconds() 

for i in range (N):
	time = (initialdelay + T[i]) / 86400. 	#The time passed from equinox till each point in orbit in days
	theta = (2*np.pi*time) / 365.256363 	#Angle between intermediate frame (s) and (epsilon) frame about common z-axis
	epsilon = 23.45 * np.pi / 180. 		#Angle between rotation axis and orbital plane normal
	x = np.cos(theta)						#components as got from document referred
	y = np.sin(theta)*np.cos(epsilon) 
	z = np.sin(theta)*np.sin(epsilon) 
	v_sun_i = np.array([x, y, z]) 	#sun vector in ECI Frame
	v_sun_i = v_sun_i/np.linalg.norm(v_sun_i.copy())
	m_si_output[i,0] = T[i] 	#first component is time
	m_si_output[i,1:4] = v_sun_i.copy();

np.savetxt("si_output.csv", m_si_output, delimiter=",") #Saves si_output to csv file

