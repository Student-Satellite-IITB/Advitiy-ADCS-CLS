import qnv
import numpy as np
import frames as fs
#Default models of sensor : (ideal sensor) measured output = actual output (no noise or bias)

def sunsensor(sat):
	v_sv_o = sat.getSun_o()
	v_q_BO = sat.getQ_BO()
	v_sv_b = qnv.quatRotate(v_q_BO,v_sv_o)
	
	return v_sv_b

def magnetometer(sat):
    v_B_o = sat.getMag_o()
    v_q_BO = sat.getQ_BO()
    v_B_b = qnv.quatRotate(v_q_BO,v_B_o)
    return v_B_b

def gps(sat):
    v_pos_m = sat.getPos() 
    v_vel_m = sat.getVel()
    time_m = sat.getTime()

    return np.hstack([v_pos_m,v_vel_m,time_m])

def gyroscope(sat):
    return sat.getW_BI_b()

def J2_propagator(sat):
    v_pos_m = sat.getPos() 
    v_vel_m = sat.getVel()
    time_m = sat.getTime()
    
    return np.hstack([v_pos_m,v_vel_m,time_m])

#Default models of controller: (no controller)

def controller(sat):
	return(np.zeros(3))

#Default models of environment: (no disturbance)
def disturbance(sat):
	return(np.zeros(3))

#Default models of estimator: (returns qBO obtained by integrator)
def estimator(sat):
    return(sat.getQ_BO())

#Default models of estimator: (returns qBO obtained by integrator)
def estimator(sat):
    return(sat.getQ_BO())
