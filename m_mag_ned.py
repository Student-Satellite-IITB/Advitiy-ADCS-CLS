import numpy as np
import datetime
from pyigrf12 import runigrf12
import constants_1U as con1U

'''
	This code takes latitude, longitude and altitude (LLA) corresponding to time data as input.
		Time: (since epoch) in seconds
		latitude: -90 to 90 degrees
		longitude: -180 to 180 degrees (-180 excluded)
		altitude: in meters
	and using standard python code for IGRF-12, output 5 column matrix
	Column 1 : time data
	       2,3,4 : Bn, Be, Bd (Magnetic field component in nED frame in nano-Tesla
	       5: 2-norm of B
'''

lla = np.genfromtxt('LLA.csv', delimiter=",") 

N = lla.shape[0]                                 
m_mag_ned = np.zeros((N,5)) #magnetic field matrix of same size as that of LLA matrix.

z1 = 0 #indicates we want magnetic field (we can also get the secular variation using 1 instead of 0 here)
z2 = 1 #indicates the height is given in km above sea level

for i in range(N):    
    print(i)
    lat = lla[i, 1]       
    lon = lla[i, 2]
    height = lla[i, 3] * 0.001 # converting altitude to km
    elapsed_t = lla[i, 0]
    e_t = datetime.timedelta(seconds = elapsed_t)

    dt = con1U.EPOCH + e_t     #present time is time of epoch + time elasped from EPOCH

    B = runigrf12(dt, z1, z2, height, lat, lon) #calling the standard function "igrf-12" which needs datetime, flag (z1 and z2) and  altitude (in km), latitude and longitude
    
    m_mag_ned[i,0]=lla[i, 0]
    m_mag_ned[i,1:5]=B         #storing returned NED magnetic field data (in nano Tesla) in matrix

np.savetxt('mag_output_ned.csv',m_mag_ned, delimiter=",")  #saving the matrix
print ("NED frame magnetic field in nano-tesla")
