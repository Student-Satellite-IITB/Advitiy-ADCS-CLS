import numpy as np
from constants_1U import RESISTANCE, INDUCTANCE, PWM_AMPLITUDE, PWM_FREQUENCY, CONTROL_STEP
import math

def leq(a, b, rel_tol=1e-12, abs_tol=1e-12):
	'''
		Function for comparison of floats, returns true if a<=b upto the tolerance
	'''
	return abs(a-b) <= abs_tol or (a < b)

def resistorPWM(v_duty_cycle,t):
	'''
		Magnetorquer is a resistor. PWM voltage is applied across it.
		Input: Duty cycle vector demanded by control law , time since start of the current signal
		Output: v_I_applied (current vector actually applied)	
	'''
	
	T = 1/PWM_FREQUENCY
	v_i_app = np.array([0.,0.,0.])
	for k in range(3):
		#condition: t < duty*T
		if math.fmod(t,T) < math.fabs(v_duty_cycle[k])*T:
			v_i_app[k] = np.sign(v_duty_cycle[k])*PWM_AMPLITUDE/RESISTANCE
		else:
			v_i_app[k] = 0.
	
	
	return v_i_app

def lrPWM(v_duty_cycle,v_i_prev,v_t_prev,t):
	'''
		Magnetorquer is a inductor-resistor in series. PWM voltage is applied across it.
		Input: duty cycle vector , electric current vector and time vector at previous edge, current time
		Output: I_applied (current actually applied)
	'''

	T = 1/PWM_FREQUENCY
	v_i_app = np.array([0.,0.,0.])
	for k in range(3):
		#condition 1: t <= duty*T 
		#condition 2: t > 0 
		if leq(math.fmod(t,T) , math.fabs(v_duty_cycle[k])*T) and (math.fmod(t,T) >0.):			
			v_i_app[k] = (PWM_AMPLITUDE - (PWM_AMPLITUDE - abs(v_i_prev[k])*RESISTANCE)*np.exp((v_t_prev[k]-t)*RESISTANCE/INDUCTANCE))
			v_i_app[k] = (np.sign(v_duty_cycle[k])/RESISTANCE)*v_i_app[k]
			
		else:
			v_i_app[k] = np.sign(v_duty_cycle[k])*v_i_prev[k]*np.exp(RESISTANCE*(v_t_prev[k]-t)/INDUCTANCE)	
	
	return v_i_app


def getCurrentList(h,v_duty_cycle):
	'''
		This functions returns current list for t=0 to t=CONTROL_STEP. Cycle is defined such that time corresponding to 
		rising or falling edge is considered in previous cycle. 
	'''
	T = 1/PWM_FREQUENCY
	N = int(CONTROL_STEP/h)

	t = np.linspace(0,CONTROL_STEP,N, endpoint=False)
	m_i_app = np.zeros((N,4))


	v_i_prev = np.array([0.,0.,0.])
	v_t_prev = np.array([0.,0.,0.])
	
	for i in range(1,N):
# =============================================================================
# 		if math.fmod(i,N/100.) == 0: 
# 			print(100*i/N ,"%"," current list done")
# 
# =============================================================================
		for k in range(3):
			#condition 1: (at high -> low)	t - duty*T > 0 and t - duty*T <=h  or
			#condition 2: (at low -> high)	0 < t and t <= h
			t1 = (math.fmod(t[i],T) - np.sign(v_duty_cycle[k])*v_duty_cycle[k]*T)  # falling edge zone
			t2 = math.fmod(t[i],T)   #rising edge zone
			
			if (t1 > 0.0 and leq(t1,h)) or (t2 > 0.0 and leq(t2,h)): 
				v_t_prev[k] = t[i-1]
				v_i_prev[k] = m_i_app[i-1,k+1]
			
		m_i_app[i,0] = t[i]
		
		m_i_app[i,1:4] = lrPWM(v_duty_cycle,v_i_prev,v_t_prev,t[i])
		
	return m_i_app