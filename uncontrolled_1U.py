import numpy as np
import satellite
import disturbance_1U as dist

from constants_1U import v_q0_BO,v_w0_BIB, R_EARTH, ALTITUDE, DELAY_STEP
from dynamics import x_dot
import frames as fs
import solver as sol
import scipy.io as sio
import os
import matplotlib.pyplot as plt
import qnv as qnv
import math



#Read position, velocity, sun-vector in ECIF from data file. temp variable for entire file
m_sgp_output_temp_i = np.genfromtxt('sgp_output_polar.csv', delimiter=",")

m_si_output_temp_i = np.genfromtxt('si_output.csv',delimiter=",")
m_light_output_temp = np.genfromtxt('light_output.csv',delimiter=",")

count = 0 # to count no. of transitions from light to eclipse
init,end = 0,0

for k in range(0,len(m_light_output_temp)-1):
	#obtain index corresponding to the start of eclipse
	l1 = m_light_output_temp[k,1]
	l2 = m_light_output_temp[k+1,1]
	if l1 ==0.5 and l2 == 0 and count == 0:	#start of eclipse
		init = k
		count = 1
		
	elif l1==0 and l2==0.5 and count == 1:	#end of eclipse
		end = k 
		
		break

#define simulation parameters

print init,end
t0 = m_sgp_output_temp_i[init,0]
tf = m_sgp_output_temp_i[end,0]	#simulation time in seconds
h = 0.1		#step size of integration in seconds
N = int((tf-t0)/MODEL_STEP)+1

#extract init to end data from temp file
m_sgp_output_i = m_sgp_output_temp_i[init:(init+N),:].copy()
m_si_output_i = m_si_output_temp_i[init:(init+N),:].copy()
m_light_output = m_light_output_temp[init:(init+N),:].copy()
print N ,'Simulation for ' ,MODEL_STEP*(N-1),'seconds'
#Initial conditions
v_q0_BO = np.array([1.,0.,0.,0.])
r=np.linalg.norm(m_sgp_output_i[0,1:4])
v_w_IO_o = np.array([0., np.sqrt(G*M_EARTH/(r)**3), 0.])
v_w0_BI_b = - v_w_IO_o.copy()

#initialize empty matrices
m_state = np.zeros((N,7))
m_q_BO = np.zeros((N,4))
m_w_BO_b = np.zeros((N,3))
m_euler_BO = np.zeros((N,3))
m_torque_dist = np.zeros((N,3))

v_q0_BI = fs.qBO2qBI(v_q0_BO,m_sgp_output_i[0,1:4],m_sgp_output_i[0,4:7])	
r=np.linalg.norm(m_sgp_output_i[0,1:4])

m_state[0,:] = np.hstack((v_q0_BI,v_w0_BI_b))
m_q_BO[0,:] = v_q0_BO.copy()
m_w_BO_b[0,:] = fs.wBIb2wBOb(v_w0_BI_b,m_q_BO[0,:],v_w_IO_o)
m_euler_BO[0,:] = qnv.quat2euler(m_q_BO[0,:])

#Make satellite object
Advitiy = satellite.Satellite(m_state[0,:],t0)
Advitiy.setControl_b(np.array([0.,0.,0.]))	#uncontrolled satellite
Advitiy.setDisturbance_b(np.array([0.,0.,0.]))

#-------------Main for loop---------------------
for  i in range(0,N-1):
	if math.fmod(i,N/100) == 0 and i>5:
		print 100*i/N 
		
	#Set satellite parameters
	
	Advitiy.setLight(m_light_output[i,1])
	Advitiy.setState(m_state[i,:])
	Advitiy.setTime(t0 + i*MODEL_STEP)
	Advitiy.setPos(m_sgp_output_i[i,1:4])
	Advitiy.setVel(m_sgp_output_i[i,4:7])
	Advitiy.setSun_i(m_si_output_i[i,1:4])
	#calculate the disturbance torques
	v_torque_gg_b = dist.ggTorqueb(Advitiy).copy()
	v_torque_aero_b = dist.aeroTorqueb(Advitiy).copy()
	v_torque_solar_b = dist.solarTorqueb(Advitiy).copy()
	
	v_torque_total_b =(v_torque_gg_b + v_torque_aero_b + v_torque_solar_b)
	Advitiy.setDisturbance_b(v_torque_total_b)
	m_torque_dist[i,:] = v_torque_total_b.copy()

	v_state_next = np.zeros((1,7))

	#Use rk4 solver to calculate the state for next step

	for j in range(0,int(MODEL_STEP/h)):		
		v_state_next = sol.rk4Quaternion(Advitiy,x_dot,h)
		Advitiy.setState(v_state_next.copy())
		Advitiy.setTime(t0 + i*MODEL_STEP + (j+1)*h)


	m_state[i+1,:] = v_state_next.copy()
	
	#Calculate observable quantities
	
	m_q_BO[i+1,:] = fs.qBI2qBO(v_state_next[0:4],m_sgp_output_i[i+1,1:4],m_sgp_output_i[i+1,4:7])
	m_w_BO_b[i+1,:] = fs.wBIb2wBOb(v_state_next[4:7],m_q_BO[i+1,:],v_w_IO_o)
	m_euler_BO[i+1,:] = qnv.quat2euler(m_q_BO[i+1,:])


#save the data files
os.chdir('Logs/')
=======
os.mkdir('polar-identity-no-dist')
os.chdir('polar-identity-no-dist')
np.savetxt('position.csv',m_sgp_output_i[:,1:4], delimiter=",")
np.savetxt('velocity.csv',m_sgp_output_i[:,4:7], delimiter=",")

np.savetxt('time.csv',m_sgp_output_i[:,0] - t0, delimiter=",")
np.savetxt('w_BOB.csv',m_w_BO_b, delimiter=",")
np.savetxt('q_BO.csv',m_q_BO, delimiter=",")
np.savetxt('state.csv',m_state, delimiter=",")
np.savetxt('euler_BO.csv',m_euler_BO, delimiter=",")
np.savetxt('disturbance.csv',m_torque_dist, delimiter=",")

print 'polar-identity-no-dist'