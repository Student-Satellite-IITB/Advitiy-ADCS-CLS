# we start code by importing modules needed in the code

import appr_actuator as act
from constants_1U import *
import controller as con
import default_blocks as defblock
import disturbance_1U as dist
from dynamics import x_dot_BO
import frames as fs
import math
import numpy as np
import os
import qnv
import satellite
import sensor
import solver as sol
from test_cases import * 

#Read position, velocity, sun-vector, light-boolean, magnetic field (in nanoTeslas) in ECIF from data file

if (orbitbool==0):
	#get data for SSO orbit
	m_sgp_output_temp_i = np.genfromtxt('sgp_output_PO.csv', delimiter=",")
	m_si_output_temp_b = np.genfromtxt('si_output_PO.csv',delimiter=",")
	m_light_output_temp = np.genfromtxt('light_output_PO.csv',delimiter=",")
	m_magnetic_field_temp_i = np.genfromtxt('mag_output_i_PO.csv',delimiter=",") 

if (orbitbool==1):
	#get data for SSO orbit
	m_sgp_output_temp_i = np.genfromtxt('sgp_output_SSO.csv', delimiter=",")
	m_si_output_temp_b = np.genfromtxt('si_output_SSO.csv',delimiter=",")
	m_light_output_temp = np.genfromtxt('light_output_SSO.csv',delimiter=",")
	m_magnetic_field_temp_i = np.genfromtxt('mag_output_i_PO.csv',delimiter=",") 

count = 0 # to count no. of transitions from light to eclipses
init,end = 0,0

for k in range(0,len(m_light_output_temp)-2): #we go to k=length-2 only because maximum index for an array is length-1 and l2 = aaray[k+1]
	#obtain index corresponding to the start of eclipse
	l1 = m_light_output_temp[k,1]
	l2 = m_light_output_temp[k+1,1]
	if l1 ==1 and l2 == 0.5 and count == 0:	#start of first eclipse
		init = k
		count = 1
		
	elif l1==0.5 and l2==1 and count == 1:	#start of second eclipse
		end = k 
		break

#define simulation parameters
print (init)
print (end)

t0 = m_sgp_output_temp_i[init,0]
tf = m_sgp_output_temp_i[end,0]	   #tf-t0 represents simulation time in seconds
h = 0.1		                       #step size of integration in seconds  
Nmodel = int((tf-t0)/MODEL_STEP)+1 #no. of time environment-cycle will run
Ncontrol = int((tf-t0)/CONTROL_STEP) #no. of time control-cycle will run

#extract init to end data from temp file
m_sgp_output_i = m_sgp_output_temp_i[init:(init+Nmodel),:].copy()
m_si_output_b = m_si_output_temp_b[init:(init+Nmodel),:].copy()
m_light_output = m_light_output_temp[init:(init+Nmodel),:].copy()
m_magnetic_field_i = m_magnetic_field_temp_i[(init-1):(init+Nmodel),:].copy()
print ((Ncontrol)*20 ,'Simulations for', Ncontrol*CONTROL_STEP, 'seconds')

#initialize empty matrices which will be needed in this simulation
v_state = np.zeros((Nmodel,7))
euler = np.zeros((Nmodel,3))
torque_dist_total = np.zeros((Nmodel,3))
torque_dist_gg = np.zeros((Nmodel,3))
torque_dist_aero = np.zeros((Nmodel,3))
torque_dist_solar = np.zeros((Nmodel,3))
torque_control = np.zeros((Nmodel,3))

#defining initial conditions
#initial state based on initial qBO and wBOB
#perfectly aligned body frame and orbit frame (v_q0_BO is initial value defined in constants)
#Body frame is not rotating wrt orbit frame (v_w0_BOB is initial value defined in constants)
v_state[0,:] = np.hstack((v_q0_BO,v_w0_BOB))                         
euler[0,:] = qnv.quat2euler(v_q0_BO)    #finding initial euler angles

#Make satellite object
Advitiy = satellite.Satellite(v_state[0,:],t0)   #t0 from line 42 of main_code

#initializing the controlTorque and measured magnetic field of satellite object as they are called in for loop before they can be set
Advitiy.setControl_b(np.array([0.,0.,0.]))		
Advitiy.setMag_b_m_c(m_magnetic_field_i[0,:]) 

#-------------Main for loop---------------------
for  i in range(0,1):  #loop for control-cycle
	
	if math.fmod(i,int(Ncontrol/100)) == 0: #we are printing percentage of cycle completed to keep track of simulation
		print (int(100*i/Ncontrol))
	#int(CONTROL_STEP/MODEL_STEP)+1
	for k in range (0,1):  #loop for environment-cycle
		#Set satellite parameters
		#state is set inside solver
		Advitiy.setPos(m_sgp_output_i[i*int(CONTROL_STEP/MODEL_STEP)+k,1:4])
		Advitiy.setVel(m_sgp_output_i[i*int(CONTROL_STEP/MODEL_STEP)+k,4:7])
		Advitiy.setLight(m_light_output[i*int(CONTROL_STEP/MODEL_STEP)+k,1])
		Advitiy.setTime(t0 + i*CONTROL_STEP + k*MODEL_STEP) #time at a cycle 

		#control
		Advitiy.setSun_i(m_si_output_b[i*int(CONTROL_STEP/MODEL_STEP)+k,1:4])
		Advitiy.setMag_i(m_magnetic_field_i[i*int(CONTROL_STEP/MODEL_STEP)+k,1:4])
		#Quest
		#magMoment_required
		#AppTorque_b
		#gyrobias

		# disturbance torque
		if (distbool == 0):
			#getting default disturbance torque (zero in our case)
			Advitiy.setDisturbance_b(defblock.disturbance(Advitiy))

		if (distbool == 1):
			#getting disturbance torque by disturbance model
			dist.ggTorqueb(Advitiy)
			dist.aeroTorqueb(Advitiy)
			dist.solarTorqueb(Advitiy)
			torque_dist_gg[i*int(CONTROL_STEP/MODEL_STEP)+k,:] = Advitiy.getggDisturbance_b()
			torque_dist_aero[i*int(CONTROL_STEP/MODEL_STEP)+k,:] = Advitiy.getaeroDisturbance_b()
			torque_dist_solar[i*int(CONTROL_STEP/MODEL_STEP)+k,:] = Advitiy.getsolarDisturbance_b()
			torque_dist_total[i*int(CONTROL_STEP/MODEL_STEP)+k,:] = torque_dist_gg[i*int(CONTROL_STEP/MODEL_STEP)+k,:] + torque_dist_aero[i*int(CONTROL_STEP/MODEL_STEP)+k,:] + torque_dist_solar[i*int(CONTROL_STEP/MODEL_STEP)+k,:]
			Advitiy.setDisturbance_b(torque_dist_total[i*int(CONTROL_STEP/MODEL_STEP)+k,:].copy())
			
		#Use rk4 solver to calculate the state for next step
		sol.rk4Quaternion(Advitiy,x_dot_BO,h)
		#storing data in matrices
		v_state[i*int(CONTROL_STEP/MODEL_STEP)+k+1,:] = Advitiy.getState()
		euler[i*int(CONTROL_STEP/MODEL_STEP)+k+1,:] = qnv.quat2euler(Advitiy.getQ_BO())

	#sensor reading
	if (sensbool == 0):
		#getting default sensor reading (zero noise in our case)
		Advitiy.setSun_b_m(defblock.sunsensor(Advitiy))
		Advitiy.setMag_b_m_p(Advitiy.getMag_b_m_c())
		Advitiy.setMag_b_m_c(defblock.magnetometer(Advitiy))
		Advitiy.setgpsData(defblock.gps(Advitiy))
		Advitiy.setOmega_m(defblock.gyroscope(Advitiy))
		Advitiy.setJ2Data(defblock.J2_propagator(Advitiy))

	if (sensbool == 1):
		#getting sensor reading from models
		Advitiy.setSun_b_m(sensor.sunsensor(Advitiy))
		Advitiy.setMag_b_m_p(Advitiy.getMag_b_m_c())
		Advitiy.setMag_b_m_c(sensor.magnetometer(Advitiy))
		Advitiy.setgpsData(sensor.GPS(Advitiy))
		Advitiy.setOmega_m(sensor.gyroscope(Advitiy))
		Advitiy.setJ2Data(sensor.J2_propagator(Advitiy))

	#Estimated quaternion
	if (estbool == 0): #qBO is same as obtained by integrator
		Advitiy.setQUEST(defblock.estimator(Advitiy))

	#if (estbool == 1): #qBO is obtained using Quest/MEKF

	#control torque
	
	if (contcons == 0):
		#getting default control torque (zero in our case)
		Advitiy.setControl_b(defblock.controller(Advitiy))
	
	#if (contcons == 1):
		#getting control torque by detumbling controller
	
	#torque applied
	
	if (actbool == 0):
		#applied torque is equal to required torque
		Advitiy.setAppTorque_b(Advitiy.getControl_b())
	
	#if (actcons == 1):
		#getting applied torque by actuator modelling (magnetic torque limitation is being considered)

#save the data files
os.chdir('Logs-Uncontrolled/')
os.mkdir('trial')
os.chdir('trial')
np.savetxt('position.csv',m_sgp_output_i[:,1:4], delimiter=",")
np.savetxt('velocity.csv',m_sgp_output_i[:,4:7], delimiter=",")
np.savetxt('time.csv',m_sgp_output_i[:,0], delimiter=",")
np.savetxt('state.csv',v_state, delimiter=",")
np.savetxt('euler.csv',euler, delimiter=",")
np.savetxt('disturbance-total.csv',torque_dist_total, delimiter=",")
np.savetxt('disturbance-gg.csv',torque_dist_gg, delimiter=",")
np.savetxt('disturbance-solar.csv',torque_dist_solar, delimiter=",")
np.savetxt('disturbance-aero.csv',torque_dist_aero, delimiter=",")
