import numpy as np
import math
import matplotlib.pyplot as plt
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
import scipy.io as sio
from pyigrf12 import runigrf12
from datetime import *
from constants_1U import *
import frames as fs


#getorbitdata.py

def filename(TT,MS,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg,fname_load):
    os.chdir('data_files/')
    #This function give variable naming to the sgp output .csv file. The filename would contain major
    #orbit elements those was used in creating the sgp output.
    
    
    try:
        a_file = open(fname_load, 'rb')
        my_input = a_file
        a_file.close
        loaded = True
        print ('File found: ' + fname_load)
        return loaded
    except IOError as e:    #if file doesnot exist already it sets flag and calculates it so you don't have to manually check
        loaded = False
        print('No existing file could be opened.')
        return loaded
    '''
    if loaded == False:
        my_input = sgp_output
         np.savetxt(fname_load,my_input, delimiter=",")
    '''

#Function getOrbitdata_TLE calls sgp_fn and create postion and velocity of satellite at different time instants.
# First It takes TLE elements (line1 and line2) from constant file then using TLE2OE it gets relevent orbital parameters  
#and then call sgp_fn. 
def getOrbitData_TLE():
    line1 = C1U.LINE1
    line2 = C1U.LINE2
    MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg,DMeanMotion,DDMeanMotion,BStar = TOE.TLE2OE(line1,line2)


    dT = C1U.dT ## dT is in seconds.
    TT = C1U.TT ## total time in minutes
    MS = C1U.MODEL_STEP ## MS in seconds
    #filename(TT,MS,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg)
    fname_load = "sgp_i_TT%g_MS%g_MMo%g_Ecc%g_Incl%g_MAnamoly%g_ArgP%g_Raan%g.csv" \
    %(TT,MS,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg)    #load from this file
    
    if filename(TT,MS,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg,fname_load) == False :
        sgp_output = np.zeros([ len(dT),7], dtype = 'float')
        sgp_output[:,0],sgp_output[:,1:4],sgp_output[:,4:7] = \
        sgp.sgp_fn(dT,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg,DMeanMotion,DDMeanMotion,BStar)
        np.savetxt(fname_load,sgp_output, delimiter=",")
        print ('csv file is created by using TLE data')


#Function getOrbitData_OrbitELement calls sgp_fn and create postion and velocity of satellite at different time instants.
# First It takes orbitalelements from cosntant file  then call sgp_fn.
def getOrbitData_OrbitElement():
    MeanMo = C1U.MeanMo 
    Eccen = C1U.Eccen 
    Incl_deg  = C1U.Incl_deg
    MeanAnamoly_deg = C1U.MeanAnamoly_deg
    ArgP =C1U.ArgP 
    RAAN_deg = C1U.RAAN_deg
    DMeanMotion = C1U.DMeanMotion
    DDMeanMotion = C1U.DDMeanMotion
    BStar = C1U.BStar
    
    dT = C1U.dT ## dT is in seconds.
    TT = C1U.TT ## total time in minutes
    MS = C1U.MODEL_STEP ## MS in seconds
    #filename(TT,MS,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg)
    fname_load = "sgp_i_TT%g_MS%g_MMo%g_Ecc%g_Incl%g_MAnamoly%g_ArgP%g_Raan%g.csv" \
    %(TT,MS,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg)    #load from this file
    

    if filename(TT,MS,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg,fname_load) == False :
        sgp_output = np.zeros([ len(dT),7], dtype = 'float')
        sgp_output[:,0],sgp_output[:,1:4],sgp_output[:,4:7] = \
        sgp.sgp_fn(dT,MeanMo,Eccen,Incl_deg,MeanAnamoly_deg,ArgP,RAAN_deg,DMeanMotion,DDMeanMotion,BStar)
        np.savetxt(fname_load,sgp_output, delimiter=",")
        print ('csv file is created by using orbital elements')

### Call getOrbitData_TLE or getOrbitData_OrbitELement
# Uncomment whichever is required
getOrbitData_TLE()
#getOrbitData_OrbitElement()



#sunmodel.py
'''
	This part of code gives sun-vector (unit vector from center of earth to center of sun) in ECI frame.
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


#lightmodel.py
'''
    This part of code generates light model output file for given orbit
    Input: orbit data file, sun vector data file
    output: N*2 array data file 
    flag = 0 eclipse
    flag = 1 light
    flag = 0.5 penumbra
    For details of model refer: 
'''


#Read SGP and sun-model data
m_sgp_output = np.genfromtxt('sgp_output.csv', delimiter=",")
m_si_output = np.genfromtxt('si_output.csv', delimiter=",")
T = m_sgp_output[:,0] #storing first element as time

N = len(T)
m_light_output = np.zeros((N,2))

r_umbra = AU*R_EARTH / (R_SUN - R_EARTH) #distance from vertex of the umbra cone to center of the earth
r_penumbra = AU*R_EARTH / (R_SUN + R_EARTH) #distance from vertex of the penumbra cone to center of the earth
alpha = np.arcsin(R_EARTH / r_umbra) #Half the aperture of the cone made by umbra (in radians)
beta = np.arcsin(R_EARTH / r_penumbra) #Half the aperture of the cone made by penumbra (in radians)

for i in range(N):
    v_pos_i = m_sgp_output[i,1:4].copy()  #position of satellite in ECI
    v_sun_i = m_si_output[i,1:4].copy()   #sun vector in ECI
    #angle between sun-vector and satellite position vector in ECI frame
    theta = np.arccos(np.dot(v_pos_i, v_sun_i) /np.linalg.norm(v_pos_i) )
    #angle between the sunvector and vector from the vertex of the umbra cone to satellite
    theta_u = np.arccos((np.dot((v_pos_i + r_umbra*v_sun_i), v_sun_i)) / np.linalg.norm(v_pos_i + r_umbra*v_sun_i))
    #angle between the negative sunvector and vector from the vertex of the penumbra cone to satellite
    theta_p =np.arccos((np.dot((v_pos_i - r_penumbra*v_sun_i), -v_sun_i)) / np.linalg.norm(v_pos_i - r_penumbra*v_sun_i))
    
    #Boolean to store whether satellite is in light or dark. 1 implies satellite is in light.
    if (theta >= np.pi/2 + alpha) & (theta_u <= alpha):
        flag = 0
    elif (theta >= np.pi/2 - beta)  & (theta_p <= beta):
        flag = 0.5
    else:
    	flag = 1

    m_light_output[i,0] = T[i]
    m_light_output[i,1] = flag
    
np.savetxt("light_output.csv", m_light_output, delimiter=',')



#getLLA.py
'''
	This part of code generates csv file for latitude, longitude, altitude (LLA) corresponding to time data.
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

#save LLA data to file
np.savetxt('LLA.csv',m_LLA, delimiter=",")
#save ecef data to the file
np.savetxt('sgp_ecef.csv',m_sgp_ecef,delimiter=",")

print ("LLA done")

#m_mag_ned.py

'''
	This part of code takes latitude, longitude and altitude (LLA) corresponding to time data as input.
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
    e_t = timedelta(seconds = elapsed_t)

    dt = EPOCH + e_t     #present time is time of epoch + time elasped from EPOCH

    B = runigrf12(dt, z1, z2, height, lat, lon) #calling the standard function "igrf-12" which needs datetime, flag (z1 and z2) and  altitude (in km), latitude and longitude
    
    m_mag_ned[i,0]=lla[i, 0]
    m_mag_ned[i,1:5]=B         #storing returned NED magnetic field data (in nano Tesla) in matrix

np.savetxt('mag_output_ned.csv',m_mag_ned, delimiter=",")  #saving the matrix
print ("NED frame magnetic field in nano-tesla")

#m_mag_eci.py

'''
	This part of code takes magnetic field in North-East-Down (NED) frame (in nT) and 
	transforms it to ECI frame.
	Output - magnetic field in ECI Frame in nanoTesla.
'''

m_mag_ned = np.genfromtxt('mag_output_ned.csv',delimiter=",")	#in nT
m_LLA = np.genfromtxt('LLA.csv',delimiter=",")	#Lat and Lon in degrees and altitude in m (check frames.latlon for details)

N = m_mag_ned.shape[0]	  #To get number of rows in array m_mag_ned
m_mag_i = np.zeros([N,4]) #no. of rows same as matrix storing NED values, 4 columns for time and 3 components of magnetic field data in ECI frame 

for k in range(N):
	T = m_mag_ned[k,0]
	m_mag_ecef = fs.ned2ecef(m_mag_ned[k,1:4].copy(),m_LLA[k,1],m_LLA[k,2]) #calling function of frame which does the conversion from NED to ECEF (Earth centered, Earth fixed frame
	m_mag_i[k,0] = T
	m_mag_i[k,1:4] = fs.ecef2ecif(m_mag_ecef.copy(),T)                      #calling function of frame which does the conversion from ECEF to ECIF 

np.savetxt('mag_output_i.csv',m_mag_i, delimiter=",")
print ("inertial magnetic field in nano-tesla")