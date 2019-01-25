import numpy as np
import unittest
import satellite
import frames 
from dynamics import x_dot_BO 
from constants_1U import m_INERTIA, v_w_IO_o
from ddt import ddt, data, unpack

@ddt
class TestDynamicsBO(unittest.TestCase):
	
	def test_zero_torque_rates_ideal_q(self):	#zero torque, zero initial wBIB, aligned frames
		'''
		For this test case set (from sixth model of Advitiy) and ALTITUDE = 700 km
		Ixx = 0.00152529
		Iyy = 0.00145111
		Izz = 0.001476
		Ixy = 0.00000437
		Iyz = - 0.00000408
		Ixz = 0.00000118
		'''
		qBO = np.array([0.0,0.,0.0,1.])
		vel = np.array([5.0,5.0,0.0])
		v_w_BI_b = np.array([0.,0.,0.])
		v_w_BO_b = frames.wBIb2wBOb(v_w_BI_b,qBO,v_w_IO_o)
		state = np.hstack((qBO,v_w_BO_b))
		mySat = satellite.Satellite(state,12.0)
		mySat.setVel(vel)
		mySat.setDisturbance_b(np.array([0.,0.,0.]))
		mySat.setControl_b(np.array([0.,0.,0.]))
		result = x_dot_BO(mySat)
		self.assertTrue(np.allclose(result,[0.,0.00053089,0.,0.,0.,0.,0.]))
	
	def test_zero_torque_ideal_q(self):		#zero torque, random initial angular velocity, aligned frames
		qBO = np.array([.0,0.,0.,1.])
		vel = np.array([5.60,-5.0,0.0])
		state = np.hstack((qBO,np.array([-0.1,0.39159385,-0.7])))
		mySat = satellite.Satellite(state,12.05)
		mySat.setVel(vel)
		mySat.setDisturbance_b(np.array([0.,0.,0.]))
		mySat.setControl_b(np.array([0.,0.,0.]))
		result = x_dot_BO(mySat)
		self.assertTrue(np.allclose(result,[-0.05,0.1957969257,-0.35, 0., 0.00453874, -0.00185081, -0.00168015]))
	
	def test_kinematics_explicitly(self):
		qBO = np.array([0.254,-0.508,np.sqrt(1-0.4**2-0.254**2-0.508**2),0.4])
		vel = np.array([5.60,-5.0,0.0])
		state = np.hstack((qBO,np.array([0.10050588,-0.05026119,-0.3014887])))
		mySat = satellite.Satellite(state,128.05)
		mySat.setVel(vel)
		mySat.setDisturbance_b(np.array([10e-10,-4e-6,-3e-5]))
		mySat.setControl_b(np.array([1e-5,1e-5,-8e-4]))
		result = x_dot_BO(mySat)
		self.assertTrue(np.allclose(result,[0.11475622,  0.06438473, -0.04115242, 0.08290271, 0.00632537, 0.00313621, -0.56264705]))
	
if __name__=='__main__':
	unittest.main(verbosity=2)