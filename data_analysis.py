import numpy as np
import matplotlib.pyplot as plt
import qnv
import os
os.chdir('Logs-Uncontrolled/PO_identity_no_dist')

time = np.genfromtxt('time.csv',delimiter=",")
v_state = np.genfromtxt('state.csv',delimiter=",")
euler = np.genfromtxt('euler.csv',delimiter=",")
#dist_total = np.genfromtxt('disturbance-total.csv',delimiter=",")
#r = np.genfromtxt('position.csv',delimiter=",")
'''
plt.plot(time,euler[:,0],label='roll')
plt.plot(time,euler[:,1],label='pitch')
plt.plot(time,euler[:,2],label='yaw')
plt.title("euler_BO in degrees")
plt.legend()
plt.show()
'''

plt.plot(time,v_state[:,0],label='q0')
plt.plot(time,v_state[:,1],label='q1')
plt.plot(time,v_state[:,2],label='q2')
plt.plot(time,v_state[:,3],label='q3')
plt.title("q_BO")
plt.legend()
plt.show()
'''
plt.plot(time,(180./np.pi)*v_w_BOB[:,0],label='wBOB_x')
plt.plot(time,(180./np.pi)*v_w_BOB[:,1],label='wBOB_y')
plt.plot(time,(180./np.pi)*v_w_BOB[:,2],label='wBOB_z')
plt.title('wBOB in degrees')
plt.legend()
plt.show()

plt.plot(time,v_q_BO[:,0],label='q0')
plt.plot(time,v_q_BO[:,1],label='q1')
plt.plot(time,v_q_BO[:,2],label='q2')
plt.plot(time,v_q_BO[:,3],label='q3')
plt.title('qBO')
plt.legend()
plt.show()

plt.plot(time,r[:,0],label='x')
plt.plot(time,r[:,1],label='y')
plt.plot(time,r[:,2],label='z')
plt.legend()
plt.title('Position in eci in m)')
plt.show()

plt.plot(time,dist_total[:,0],label="t_x")
plt.plot(time,dist_total[:,1],label="t_y")
plt.plot(time,dist_total[:,2],label="t_z")
plt.legend()
plt.title('disturbance torque')
plt.show()

plt.plot(time,(180./np.pi)*v_state[:,4],label='wBIBx')
plt.plot(time,(180./np.pi)*v_state[:,5],label='wBIBy')
plt.plot(time,(180./np.pi)*v_state[:,6],label='wBIBz')
plt.title("wBIB in degrees")
plt.legend()
plt.show()
'''
