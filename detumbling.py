import numpy as np
import satellite
import disturbance_1U as dist
from constants_1U import v_q0_BO, MODEL_STEP, LINE1, LINE2, G, m_INERTIA, M_EARTH, No_Turns, v_A_Torquer, PWM_AMPLITUDE, PWM_FREQUENCY, RESISTANCE
from dynamics import x_dot_BI
import frames as fs
import solver as sol
import os
import qnv
import math
import detumbling_con as detcon
#import actuator
from TorqueApplied import ctrlTorqueToVoltage,currentToTorque,I
import actuator as act

#Read position, velocity, sun-vector, magnetic field (in nanoTeslas) in ECIF from data file

m_sgp_output_temp_i = np.genfromtxt('sgp_output.csv', delimiter=",")
m_si_output_temp_b = np.genfromtxt('si_output.csv',delimiter=",")
m_light_output_temp = np.genfromtxt('light_output.csv',delimiter=",")
m_magnetic_field_temp_i = np.genfromtxt('mag_output_i.csv',delimiter=",") 

count = 0 # to count no. of transitions from light to eclipses
init,end = 0,0

for k in range(0,len(m_light_output_temp)-1):
	#obtain index corresponding to the start of eclipse
	l1 = m_light_output_temp[k,1]
	l2 = m_light_output_temp[k+1,1]
	if l1 ==0.5 and l2 == 0 and count == 0:	#start of first eclipse
		init = k
		count = 1
		
	elif l1==0.5 and l2==0 and count == 1:	#start of second eclipse
		end = k 
		break

#define simulation parameters
print (init)
print (end)

t0 = m_sgp_output_temp_i[init,0]
tf = m_sgp_output_temp_i[end,0]	#simulation time in seconds
#h = 0.1		#step size of integration in seconds  #declared in constants_1U.py
N = int((tf-t0)/MODEL_STEP)+1

#extract init to end data from temp file
m_sgp_output_i = m_sgp_output_temp_i[init:(init+N),:].copy()
m_si_output_b = m_si_output_temp_b[init:(init+N),:].copy()
m_light_output = m_light_output_temp[init:(init+N),:].copy()
m_magnetic_field_i = m_magnetic_field_temp_i[(init-1):(init+N),:].copy()
print (N ,'Simulation for ' ,MODEL_STEP*(N-1),'seconds')

#initialize empty matrices
v_state = np.zeros((N,7))
v_q_BO = np.zeros((N,4))
v_w_BOB = np.zeros((N,3))
euler = np.zeros((N,3))
torque_dist = np.zeros((N,3))

v_current_p = np.zeros(3)
v_current_req = np.zeros(3)
v_q0_BI = fs.qBO_2_qBI(v_q0_BO,m_sgp_output_i[0,1:4],m_sgp_output_i[0,4:7])	
r=np.linalg.norm(m_sgp_output_i[0,1:4])
v_w0_BIB = -np.array([0., np.sqrt(G*M_EARTH/(r)**3), 0.])
v_state[0,:] = np.hstack((v_q0_BI,v_w0_BIB))
v_q_BO[0,:] = v_q0_BO
v_w_BOB[0,:] = fs.wBIb2wBOb(v_w0_BIB,v_q_BO[0,:],(-v_w0_BIB))
euler[0,:] = qnv.quat2euler(v_q_BO[0,:])

#Make satellite object
Advitiy = satellite.Satellite(v_state[0,:],t0)
Advitiy.setControl_b(np.array([0.,0.,0.]))	
Advitiy.setDisturbance_i(np.array([0.,0.,0.]))

#-------------Main for loop---------------------
for  i in range(0,N-1):
	
	if math.fmod(i,int(N/100)) == 0:
		print (int(100*i/N)) 
	
	#Set satellite parameters
	Advitiy.setLight(m_light_output[i,1])
	Advitiy.setState(v_state[i,:])
	Advitiy.setTime(t0 + i*MODEL_STEP)
	Advitiy.setPos(m_sgp_output_i[i,1:4])
	Advitiy.setVel(m_sgp_output_i[i,4:7])
	Advitiy.setSun_i(m_si_output_b[i,1:4])
	#v_magnetic_field_b_p=sen
	# obtain these data from magmeter modelling
	#print (m_magnetic_field_i[i,1:4])
	if i==0:
		v_magnetic_field_b_p=qnv.quatRotate(v_state[0,0:4],m_magnetic_field_i[i,1:4])
	else:
		v_magnetic_field_b_p=qnv.quatRotate(v_state[i-1,0:4],m_magnetic_field_i[i,1:4])
	v_magnetic_field_b_c=qnv.quatRotate(Advitiy.getQ(),m_magnetic_field_i[i+1,1:4])
	Advitiy.setMag_b_m_p(v_magnetic_field_b_p) 		
	Advitiy.setMag_b_m_c(v_magnetic_field_b_c)
	Advitiy.setMag_i(m_magnetic_field_i[i+1,1:4])

	v_torque_gg_b = dist.ggTorqueb(Advitiy).copy()
	v_torque_aero_b = dist.aeroTorqueb(Advitiy).copy()
	v_torque_solar_b = dist.solarTorqueb(Advitiy).copy()
	v_torque_total_b =(v_torque_gg_b + v_torque_aero_b + v_torque_solar_b)
	Advitiy.setDisturbance_i(v_torque_total_b)
	torque_dist[i,:] = v_torque_total_b.copy()
	
	if math.fmod(i,20) == 0:
		v_magMoment = detcon.magMoment(Advitiy)
	v_magnetic_field_b = qnv.quatRotate(Advitiy.getQ(),Advitiy.getMag_i()) * 10**(-9)
	v_control_torque_b = np.cross(v_magMoment,v_magnetic_field_b)
	Advitiy.setControl_b(v_control_torque_b)

	v_state_next = np.zeros((1,7))

    if i%20==0:
        voltage=ctrlTorqueToVoltage(Advitiy)
        v_duty_cycle=voltage/PWM_AMPLITUDE
        m_current_list = act.getCurrentList(h,v_duty_cycle)  #for getting  PWM current list for a CONTROL_STEP
# =============================================================================
#         m_current_list=I(voltage)    # for getting DC current list for a CONTROL_STEP
# =============================================================================
        v_app_torque_b=currentToTorque(m_current_list,Advitiy)
        for k in range(0,v_app_torque_b.shape[0]):
            Advitiy.setAppTorque_b(v_app_torque_b[k].copy())
            Advitiy.setTime(t0 + i*MODEL_STEP + (k+1)*h)
            
	#Use rk4 solver to calculate the state for next step
	for j in range(0,int(MODEL_STEP/h)):		
		v_state_next = sol.rk4Quaternion(Advitiy,x_dot,h)
		Advitiy.setState(v_state_next.copy())
		Advitiy.setTime(t0 + i*MODEL_STEP + (j+1)*h)

	v_state[i+1,:] = v_state_next.copy()
	
	#Calculate observable quantities
	v_q_BO[i+1,:] = fs.qBI2qBO(v_state_next[0:4],m_sgp_output_i[i+1,1:4],m_sgp_output_i[i+1,4:7])
	r=np.linalg.norm(m_sgp_output_i[i+1,1:4])
	v_w0_BIB = -np.array([0., np.sqrt(G*M_EARTH/(r)**3), 0.])
	v_w_BOB[i+1,:] = fs.wBIb2wBOb(v_state_next[4:7],v_q_BO[i+1,:],(-v_w0_BIB))
	euler[i+1,:] = qnv.quat2euler(v_q_BO[i+1,:])

#save the data files
os.chdir('Logs-Detumbling/')
os.mkdir('trial')
os.chdir('trial')
np.savetxt('position.csv',m_sgp_output_i[init:end+1,1:4], delimiter=",")
np.savetxt('velocity.csv',m_sgp_output_i[init:end+1,4:7], delimiter=",")
np.savetxt('time.csv',m_sgp_output_i[init:end+1,0] - t0, delimiter=",")
np.savetxt('w_BOB.csv',v_w_BOB, delimiter=",")
np.savetxt('q_BO.csv',v_q_BO, delimiter=",")
np.savetxt('state.csv',v_state, delimiter=",")
np.savetxt('euler.csv',euler, delimiter=",")
#np.savetxt('disturbance.csv',torque_dist, delimiter=",")
np.savetxt('Moment of Inertia', m_INERTIA, delimiter=",")
