#Code for transformation of vector from one reference frame to another

import numpy as np
from constants_1U import W_EARTH, EPOCH, EQUINOX, STEPRUT
from math import radians, sin, cos, acos, pi
import qnv as qnv

def latlon(v_x):
	#get the latitude and longitude in degrees given position in ECEF 
	#Input: x is position in ecef frame
	#Output: latitude and longitude in degres (tuple)
	#latitude ranges [0,90] in north hemisphere and [0,-90] in south hemisphere

	lat = sgn(v_x[2])*(acos(((v_x[0]**2 + v_x[1]**2)**0.5)/((v_x[0]**2 + v_x[1]**2 + v_x[2]**2)**0.5)))*90./(pi/2.) 
	   
	# longitude calculation given position, lon is longitude
	#ranges from (-pi,pi]
	if v_x[1]==0:
		if v_x[0]>=0 :
			lon = 0.0
		else:
			lon = 180.0

	else:
		lon = sgn(v_x[1])*acos(v_x[0]/((v_x[0]**2 + v_x[1]**2)**0.5))*90./(pi/2)      

	#x axis is intersection of 0 longitude and 0 latitude

	return lat,lon

def sgn(x):
	#signum function
	if(x==0):
		y = 0
	elif (x>0):
		y = 1
	else:
		y = -1

	return y

def ecif2ecef(v_x_i,t):
	#Input ecif cector and time since epoch in seconds
	# output ecef vector
	ut_sec = (EPOCH - EQUINOX).total_seconds() + t # universal time vector in sec 
	theta = W_EARTH*ut_sec #in radian
	m_DCM = np.array([[cos(theta), sin(theta), 0.], [-1*sin(theta), cos(theta),0.], [0.,0.,1.]])
	v_x_e = np.dot(m_DCM,v_x_i)
	
	return v_x_e

def ecif2ecefR(t):

	ut_sec = (EPOCH - EQUINOX).total_seconds() + t # universal time vector in sec
	st_sec = STEPRUT*ut_sec   #sidereal time vector in sec

	phi = st_sec*W_EARTH           # sidereal time vector in rad

	TEI = np.array([[cos(phi),sin(phi),0.],[-sin(phi),cos(phi),0.],[0.,0.,1.]])

	return TEI


def ecef2ecif(v_x_e,t):
	#Input ecef cector and time since epoch in seconds
	# output ecif vector
	ut_sec = (EPOCH - EQUINOX).total_seconds() + t # universal time vector in sec
	theta = W_EARTH*ut_sec #in radian
	m_DCM = np.array([[cos(theta), -1*sin(theta), 0.], [sin(theta), cos(theta),0.],[ 0.,0.,1.]])
	v_x_i = np.dot(m_DCM,v_x_e)
	
	return v_x_i

def ecif2orbit(v_pos_i,v_vel_i,v_x_i):
	#Input: v_pos_i is position in eci frame , v_vel_i is velocity in eci frame, v_x_i is vector to be transformed
	#output: vector components in orbit frame
	z = -v_pos_i/np.linalg.norm(v_pos_i)
	y = np.cross(v_vel_i,v_pos_i)/np.linalg.norm(np.cross(v_vel_i,v_pos_i))
	x = np.cross(y,z)/np.linalg.norm(np.cross(y,z))
	m_DCM_OI = np.array([x,y,z])
	v_x_o = np.dot(m_DCM_OI,v_x_i)

	return v_x_o

def qBI2qBO(v_q_BI,v_pos_i,v_vel_i):
	#input: Quaternion corresponding to transformation of components of vector in ecif to components in body frame, 
	#		position and velocity in ecif
	#output: Quaternion corresponding to transformation of components of vector in orbit frame to components in body frame
	#Orbit frame def:	z_o = nadir(opposite to satellite position vector) y_o: cross(v,r) x_o: cross(y,z) 

	z = -v_pos_i/np.linalg.norm(v_pos_i)
	y = np.cross(v_vel_i,v_pos_i)/np.linalg.norm(np.cross(v_vel_i,v_pos_i))
	x = np.cross(y,z)/np.linalg.norm(np.cross(y,z))
	m_DCM_OI = np.array([x,y,z])
	v_q_IO = qnv.rotm2quat(np.transpose(m_DCM_OI))
	v_q_BO = qnv.quatMultiplyNorm(v_q_BI,v_q_IO)
	if v_q_BO[3] < 0.:
		v_q_BO = -1.*v_q_BO.copy()
	v_q_BO = v_q_BO/np.linalg.norm(v_q_BO.copy())	

	return v_q_BO

def qBO2qBI(v_q_BO,v_pos_i,v_vel_i):
	#input: Quaternion corresponding to transformation of components of vector in orbit frame to components in body frame, 
	#		position and velocity in ecif
	#output: Quaternion corresponding to transformation of components of vector in ecif to components in body frame
	#Orbit frame def:	z_o = nadir(opposite to satellite position vector) y_o: cross(v,r) x_o: cross(y,z)
	z = -v_pos_i/np.linalg.norm(v_pos_i)
	y = np.cross(v_vel_i,v_pos_i)/np.linalg.norm(np.cross(v_vel_i,v_pos_i))
	x = np.cross(y,z)/np.linalg.norm(np.cross(y,z))
	m_DCM_OI = np.array([x,y,z])

	v_q_OI = qnv.rotm2quat(m_DCM_OI)
	v_q_BI = qnv.quatMultiplyNorm(v_q_BO,v_q_OI)
	if v_q_BI[3] < 0. :
		v_q_BI = -1.*v_q_BI.copy()
	v_qBI = v_q_BI/np.linalg.norm(v_q_BI.copy())	
	return v_q_BI

def wBIb2wBOb(v_w_BI_b,v_q_BO,v_w_IO_o):
	#input: angular velocity of body wrt ecif in body frame, unit quaternion which rotates orbit vector to body frame
	#		angular velocity of ecif wrt orbit frame in orbit frame
	#output: angular velocity of body frame wrt orbit frame in body frame
	
	return v_w_BI_b + qnv.quatRotate(v_q_BO,v_w_IO_o)

def ned2ecef(v,lat,lon):
	#rotate vector from North-East-Down frame to ecef
	#Input:
	#	lat: latittude in degrees ranges from -90 to 90
	#	lon: longitude in degrees ranges from (-180,180]	 
	if lat == 90 or lat == -90:
		raise ValueError('Latittude value +/-90 occured. NED frame is not defined at north and south pole !!')
	theta = -lat + 90. #in degree, polar angle (co-latitude)

	if lon<0:
		phi = 360. - abs(lon)
	else:
		phi = lon #in degree, azimuthal angle
	theta = radians(theta)
	phi = radians(phi)

	m_DCM_n2e = np.array([[ -cos(theta)*cos(phi),	-sin(phi),	-sin(theta)*cos(phi)],
						[	-cos(theta)*sin(phi),	cos(phi),	-sin(theta)*sin(phi)],
						[	sin(theta),	0.0,	-cos(theta)]])
	
	y = np.dot(m_DCM_n2e,v)

	return y

def wBOb2wBIb(v_w_BO_b,v_q_BO,v_w_IO_o):
	#input: angular velocity of body wrt orbit in body frame, unit quaternion which rotates orbit vector to body frame
	#		angular velocity of ecif wrt orbit frame in orbit frame
	#output: angular velocity of body frame wrt eci frame in body frame
	
	return v_w_BO_b - qnv.quatRotate(v_q_BO,v_w_IO_o)

