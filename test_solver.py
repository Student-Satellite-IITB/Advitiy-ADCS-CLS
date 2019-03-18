from solver import updateStateTimeRK4
import numpy as np
import frames
import unittest
import satellite
from dynamics import x_dot_BO
from constants_1U import G, M_EARTH, v_w_IO_o
from ddt import data, ddt

def f_constant(sat):
	cons = np.array((1.,1.,1.,1.,1.,1.,1.))
	return cons
@ddt
class TestSolver(unittest.TestCase):
	
	def test_constant(self):
		t0 = 0.0
		mySat = satellite.Satellite(np.array([1.0,0.,0.,0.,0.,0.,0.]),t0)
		t = 10
		mySat.setPos(np.array([1e6,0.,0.]))
		mySat.setVel(np.array([5.0,5.0,0.0]))
		h = 0.1
		for i in range(0,int(t/h)):
			updateStateTimeRK4(mySat,f_constant,h)
			mySat.setTime(t0+(i+1)*h)
		self.assertTrue(np.allclose([10.0,10.0,10.0],mySat.getW_BO_b()))
		
	@data(0.001,0.005,0.05,0.1,0.5,1,5,10)	
	def test_dynamics(self,value):
		t0 = 0.
		h = value
		
		v_q_BO = np.array([0.4,0.254,-0.508,0.71931912])
		
		v_w_BO_b = frames.wBIb2wBOb(np.array([0.1,-0.05,-0.3]),v_q_BO,v_w_IO_o)
		
		mySat1 = satellite.Satellite(np.hstack((v_q_BO, v_w_BO_b)),t0)
		
		mySat1.setPos(np.array([1e6,53e5,0.]))
		mySat1.setVel(np.array([5.60,-5.0,0.0]))
		mySat1.setDisturbance_b(np.array([10e-10,-4e-6,-3e-5]))
		mySat1.setControl_b(np.array([1e-5,1e-5,-8e-4]))
		
		updateStateTimeRK4(mySat1,x_dot_BO,h)
		x1 = mySat1.getState()

		mySat2 = satellite.Satellite(np.hstack((v_q_BO, v_w_BO_b)),t0)
		
		mySat2.setPos(np.array([1e6,53e5,0.]))
		mySat2.setVel(np.array([5.60,-5.0,0.0]))
		mySat2.setDisturbance_b(np.array([10e-10,-4e-6,-3e-5]))
		mySat2.setControl_b(np.array([1e-5,1e-5,-8e-4]))
		
		mySat2.setState(np.hstack((v_q_BO,v_w_BO_b)))		
		k1 = h*x_dot_BO(mySat2)
		mySat2.setState(np.hstack((v_q_BO,v_w_BO_b))+0.5*k1)
		k2 = h*x_dot_BO(mySat2)
		mySat2.setState(np.hstack((v_q_BO,v_w_BO_b))+0.5*(k2))
		k3 = h*x_dot_BO(mySat2)
		mySat2.setState(np.hstack((v_q_BO,v_w_BO_b))+k3)
		k4 = h*x_dot_BO(mySat2)
		
		error_state = np.hstack((v_q_BO,v_w_BO_b)) + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)
		if (error_state[3]<0):
			error_state[0:4]=-1*error_state[0:4]
		error_state[0:4]=error_state[0:4]/np.linalg.norm(error_state[0:4])
		self.assertTrue(np.allclose(error_state,x1))
	
if __name__=="__main__":
	unittest.main(verbosity=2)