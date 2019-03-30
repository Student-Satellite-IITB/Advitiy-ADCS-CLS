import numpy as np
import math
import frames as fs

def quatInv(v_q1): #to get inverse of a quaternion
	v_qi = np.zeros(4)
	v_qi[3] = v_q1[3]
	v_qi[0:3] = -1*v_q1[0:3].copy()
	return v_qi

def quatMultiplyNorm(v_q1,v_q2): #returns quaternion product (product is a unit quaternion)

	a1 = v_q1[3].copy()
	a2 = v_q2[3].copy()
	
	v_b1 = (v_q1[0:3].copy())
	v_b2 = (v_q2[0:3].copy())
	
	a = a1*a2 - np.dot(v_b1,v_b2)
	v_b = a1*v_b2 + a2*v_b1 - np.cross(v_b1,v_b2)
	v_q = np.hstack((v_b,a))
	v_q = v_q/np.linalg.norm(v_q)
	return v_q

def quatMultiplyUnnorm(v_q1,v_q2): #returns quaternion product (product is not a unit quaternion)

	a1 = v_q1[3].copy()
	a2 = v_q2[3].copy()
	
	v_b1 = (v_q1[0:3].copy())
	v_b2 = (v_q2[0:3].copy())
	
	a = a1*a2 - np.dot(v_b1,v_b2)
	v_b = a1*v_b2 + a2*v_b1 - np.cross(v_b1,v_b2)

	v_q = np.hstack((v_b,a))
	return v_q

def quatRotate(v_q,v_x): 
	#Transforms components of vector x in frame A to components in frame B if the input quaternion is q_BA	
	if np.count_nonzero(v_x) == 0:
		return v_x
	v_qi = quatInv(v_q)
	v_y = np.hstack((v_x.copy(),[0.]))
	v_y = quatMultiplyUnnorm(v_q,v_y)
	v_y = quatMultiplyUnnorm(v_y,v_qi)
	v_x2 = v_y[0:3]
	return v_x2

def quatDerBI(v_q,v_w): 	# w is angular velocity of body wrt inertial frame in body frame. 

						#q transforms inertial frame vector to body frame

	m_W = np.array([[0.,v_w[2],-v_w[1],v_w[0]],[-v_w[2],0.,v_w[0],v_w[1]],[v_w[1],-v_w[0],0.,v_w[2]],[-v_w[0],-v_w[1],-v_w[2],0.]])
	v_q_dot = 0.5*np.dot(m_W,v_q)

	return v_q_dot

def quatDerBO(v_q,v_w):   #q transforms orbit frame vector to body frame 
						  # w is angular velocity of body wrt orbit frame in body frame)

	m_W = np.array([[0.,v_w[2],-v_w[1],v_w[0]],[-v_w[2],0.,v_w[0],v_w[1]],[v_w[1],-v_w[0],0.,v_w[2]],[-v_w[0],-v_w[1],-v_w[2],0.]])
	v_q_dot = 0.5*np.dot(m_W,v_q)
	return v_q_dot

def rotm2quat(m_A): #returns a quaternion whose scalar part is positive to keep angle between -180 to +180 deg.
					#formula from pg 97 of textbook by Junkins
	q4 = 1 + np.trace(m_A)
	q1 = 1 + m_A[0,0] - m_A[1,1] - m_A[2,2]
	q2 = 1 - m_A[0,0] + m_A[1,1] - m_A[2,2]
	q3 = 1 - m_A[0,0] - m_A[1,1] + m_A[2,2]
	qm = max(q1,q2,q3,q4)
	if(qm==q4):
		q4 = math.sqrt(q4)/2
		q1 = (m_A[1,2] - m_A[2,1])/(4*q4)
		q2 = (m_A[2,0] - m_A[0,2])/(4*q4)
		q3 = (m_A[0,1] - m_A[1,0])/(4*q4)

	elif(qm==q1):
		q1 = math.sqrt(q1)/2
		q4 = (m_A[1,2] - m_A[2,1])/(4*q1)
		q2 = (m_A[0,1] + m_A[1,0])/(4*q1)
		q3 = (m_A[0,2] + m_A[2,0])/(4*q1)

	elif(qm==q2):
		q2 = math.sqrt(q2)/2
		q4 = (m_A[2,0] - m_A[0,2])/(4*q2)
		q1 = (m_A[0,1] + m_A[1,0])/(4*q2)
		q3 = (m_A[1,2] + m_A[2,1])/(4*q2)

	else: 
		q3 = math.sqrt(q3)/2
		q4 = (m_A[0,1] - m_A[1,0])/(4*q3)
		q1 = (m_A[0,2] + m_A[2,0])/(4*q3)
		q2 = (m_A[1,2] + m_A[2,1])/(4*q3)

	v_q = np.array([q1,q2,q3,q4])
	v_q = v_q/np.linalg.norm(v_q)
	if (fs.sgn(v_q[3])==-1):
		v_q = -1 * v_q
	return v_q

def quat2rotm(v_q): #given a quaternion it returns a rotation matrix
	q1 = v_q[0]  
	q2 = v_q[1]
	q3 = v_q[2]
	q4 = v_q[3]

	m_M1 = 2* np.array([[-q2**2 - q3**2,q1*q2,q1*q3],[q1*q2,-q1**2 - q3**2,q2*q3],[q1*q3,q2*q3,-q1**2-q2**2]])
	m_M2 = -2*q4*np.array([[0,-q3,q2],[q3,0,-q1],[-q2,q1,0]])
	m_M3 = np.identity(3)  #norm(v_q)=1 
	return m_M1 + m_M2 + m_M3


def quat2euler(v_q):
	#input quaternion
	#output euler angles: roll, pitch, yaw in degrees
	m_M = quat2rotm(v_q)
	yaw = math.atan2(m_M[0,1],m_M[0,0])
	pitch = -math.asin(m_M[0,2])
	roll = math.atan2(m_M[1,2],m_M[2,2])
	return (180./np.pi)*np.array([roll,pitch,yaw])
