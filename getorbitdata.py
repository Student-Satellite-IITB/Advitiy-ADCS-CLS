import numpy as np
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from datetime import *
from constants_1U import LINE1, LINE2, EPOCH

satellite = twoline2rv(LINE1, LINE2, wgs72) #wgs72 is a particular model used by sgp4

delay = 0.1 #Enter the delay between each consecutive data points in seconds
totaltime = 100.*60. #Enter the time for which SGP Data should be created from the initial time 
N = int(totaltime/delay) #Number of iterations for which SGP4 module will be called
delay = timedelta(seconds=delay) #Converts the delay into a form that can be added to datetime objects
sgp_output = np.zeros([N,7])
time = EPOCH
for i in range (N):
	time = time+delay
	sgp_output[i,0] = i*0.1 	#Stores the time after launchdate for which positon and velocity is calculated
	sgp_output[i,1:4]=satellite.propagate(time.year,time.month,time.day,time.hour,time.minute,time.second)[0] #Stores the position
	sgp_output[i,4:7]=satellite.propagate(time.year,time.month,time.day,time.hour,time.minute,time.second)[1] #Stores the velocity, Reshape is necessary to match dimensions
np.savetxt("sgp_output.csv", sgp_output, delimiter=",") #Saves sgp_output to csv file
print type(sgp_output)
print np.shape(sgp_output)

print sgp_output[0,:]
print sgp_output[3600,:]
print sgp_output[7200,:]
print sgp_output[10800,:]
print sgp_output[14400,:]