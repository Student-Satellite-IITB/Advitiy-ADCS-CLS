import numpy as np
import frames as fs
import disturbance_1U as dist
import unittest
from satellite import Satellite
from constants_1U import m_INERTIA,v_w_IO_o
from ddt import data, ddt , unpack

@ddt
class TestGravityGradient(unittest.TestCase):
	#Test the data type of output
	def test_gg_data_type(self):
		qBI = np.array([0.,0.,0.,1.])
		wBIb = np.array([0.,0.,0.])
		v_pos_i = np.array([7070e3,0.,0.])
		v_vel_i = np.array([2.0e3,2.8,-73.2])
		qBO = fs.qBI2qBO(qBI,v_pos_i,v_vel_i)
		wBOB = fs.wBIb2wBOb(wBIb,qBO,v_w_IO_o)
		sat = Satellite(np.hstack((qBO,wBOB)),13.)
		sat.setPos(np.array([7070e3,0.,0.]))  
		sat.setVel(v_vel_i)
		
		dist.ggTorqueb(sat)
		result = sat.getggDisturbance_b()        
		self.assertEqual(type(result),np.ndarray)
	
	#Take position vector to be eigen vector of m_INERTIA along with q =[1,0,0,0] thus torque will be zero
	l, w = np.linalg.eig(m_INERTIA)
	@data(7070e3*w[:,0],7070e3*w[:,1],7070e3*w[:,2])
	def test_inertia_eigenvec(self,value):
		
		qBI = np.array([0.,0.,0.,1.])
		wBIb = np.array([0.1,-0.02,-0.2])
		v_pos_i = value
		v_vel_i = np.array([5.60,-5.0,0.0])
		qBO = fs.qBI2qBO(qBI,v_pos_i,v_vel_i)
		wBOB = fs.wBIb2wBOb(wBIb,qBO,v_w_IO_o)
		sat = Satellite(np.hstack((qBO,wBOB)),13.)
		sat.setPos(value)
		sat.setVel(v_vel_i)
		dist.ggTorqueb(sat)
		result = sat.getggDisturbance_b()
		self.assertTrue(np.allclose(result,[0.,0.,0.]))      
	
class TestAeroDrag(unittest.TestCase):
	#test the data type of the output
	
	def test_aero_type(self):
		qBI = np.array([-np.sqrt(0.5),0.,0.,np.sqrt(0.5)])
		wBIb = np.array([0.1,0.23,0.])
		v_pos_i = np.array([0.,0.,7e6])
		v_vel_i = np.array([0,2e3,6e3])
		qBO = fs.qBI2qBO(qBI,v_pos_i,v_vel_i)
		wBOB = fs.wBIb2wBOb(wBIb,qBO,v_w_IO_o)
		sat = Satellite(np.hstack((qBO,wBOB)),13.)
		sat.setPos(np.array([0.,0.,7e6]))		
		sat.setVel(np.array([0,2e3,6e3]))
		dist.aeroTorqueb(sat)
		result = sat.getaeroDisturbance_b()
			   
		self.assertEqual(type(result),np.ndarray)

	def test_aero_value(self):
		qBI = np.array([-np.sqrt(0.5),0.,0.,np.sqrt(0.5)])
		wBIb = np.array([0.1,0.23,0.])
		v_pos_i = np.array([0.,0.,7e6])
		v_vel_i = np.array([0,2e3,6e3])
		qBO = fs.qBI2qBO(qBI,v_pos_i,v_vel_i)
		wBOB = fs.wBIb2wBOb(wBIb,qBO,v_w_IO_o)
		sat = Satellite(np.hstack((qBO,wBOB)),13.)
		sat.setPos(np.array([0.,0.,7e6]))
		sat.setVel(np.array([0,2e3,6e3]))
		dist.aeroTorqueb(sat)
		result = sat.getaeroDisturbance_b()
		self.assertTrue(np.allclose(result, [2.99654080e-10,-2.57065600e-11,-7.71196800e-11]))
	   
		
class TestSolarTorque(unittest.TestCase):
		#test the data type of the output
	def test_solar_torque_type(self):
		qBI = np.array([0.,0.,0.,1.])
		wBIb = np.array([0.,0.,0.])
		v_pos_i = np.array([7070e3,0.,0.])
		v_vel_i = np.array([0,2e3,6e3])
		qBO = fs.qBI2qBO(qBI,v_pos_i,v_vel_i)
		wBOB = fs.wBIb2wBOb(wBIb,qBO,v_w_IO_o)
		sat = Satellite(np.hstack((qBO,wBOB)),13.)
		sat.setPos(np.array([7070e3,0.,0.]))
		v_sv_i=np.array([1.0,0.0,0.0])
		sat.setSun_i(v_sv_i)
		sat.setLight(1)  
		sat.setPos(v_pos_i)
		sat.setVel(v_vel_i)
		dist.solarTorqueb(sat)
		result = sat.getsolarDisturbance_b()        
		self.assertEqual(type(result),np.ndarray)
	
	def test_solar_torque_value(self):
		qBI = np.array([0.,0.,0.,1.])
		wBIb = np.array([0.1,-0.02,-0.2])
		v_pos_i = np.array([7070e3,0.,0.])
		v_vel_i = np.array([0,2e3,6e3])
		qBO = fs.qBI2qBO(qBI,v_pos_i,v_vel_i)
		wBOB = fs.wBIb2wBOb(wBIb,qBO,v_w_IO_o)
		sat = Satellite(np.hstack((qBO,wBOB)),13.)
		v_sv_i=np.array([1.0,0.0,0.0])            #sun vector in eci frame
		sat.setSun_i(v_sv_i)
		sat.setLight(1)
		sat.setPos(v_pos_i)
		sat.setVel(v_vel_i)
		dist.solarTorqueb(sat)
		result = sat.getsolarDisturbance_b()        
		self.assertTrue(np.allclose(result,[  0.00000000e+00,-3.66624000e-11,3.17376000e-10]))

if __name__=="__main__":
	unittest.main(verbosity=2)