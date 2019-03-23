import frames as fr 	#module to be tested
import unittest	#testing library
import numpy as np
from ddt import ddt,file_data,unpack,data
from math import sin,cos
from constants_1U import W_EARTH, EPOCH, EQUINOX
import qnv as qnv

@ddt
class TestLatLon(unittest.TestCase):
	@file_data("test-data/test_latlon.json")	#file for data sets
	@unpack
	def test_cases_file(self,value):
		given,expected = np.asarray(value[0]),value[1]	#asarray to convert list to array
		result = fr.latlon(given)
		self.assertTrue(np.allclose([result[0],result[1]],expected))	
		self.assertEqual(type(result),tuple)

	@data(1,2,3,4,5,6)
	def test_z_symmetry_lat(self,value):
		#self.assertEqual(2,2)
		x1 = np.array([sin(value),cos(value),13.])
		self.assertAlmostEqual(fr.latlon(x1)[0],fr.latlon(np.array([1.,0.,13.]))[0])

	@data(7e3,7.1e3,7.2e3,7.3e3,7.4e3)
	def test_independnt_of_radius(self,radius):
		x1 = np.array([1.,2.,-3.])
		self.assertAlmostEqual(fr.latlon(radius*x1)[0],fr.latlon(x1)[0])
		self.assertAlmostEqual(fr.latlon(radius*x1)[1],fr.latlon(x1)[1])

@ddt
class TestEcif2Ecef(unittest.TestCase):
	@data(100.0,210.,323.0,452.0,580.90)
	def test_orthogonality_of_matrix(self,value):	#Test to check that the rotation matrix is orthogonal
		t = (EPOCH - EQUINOX).total_seconds() + value

		v_x = fr.ecif2ecef(np.array([1.0,0.0,0.0]),t)
		v_y = fr.ecif2ecef(np.array([0.0,1.0,0.0]),t)
		v_z = fr.ecif2ecef(np.array([0.0,0.0,1.0]),t)
		self.assertEqual(type(v_x),np.ndarray)	

		m = np.array([v_x,v_y,v_z])
		mT = m.transpose()
		I = np.dot(m,mT)
		self.assertTrue(np.allclose(I[0,:],[1.,0.,0.]))
		self.assertTrue(np.allclose(I[1,:],[0.,1.,0.]))
		self.assertTrue(np.allclose(I[2,:],[0.,0.,1.]))

	def test_at_zero_time(self):	#Test when ecef and ecif are aligned
		t = (EQUINOX - EPOCH).total_seconds()	
		u = np.array([5.0e3,-2.30e3,-5.89e3])
		v = fr.ecif2ecef(u,t)
		self.assertTrue(np.allclose(v,u))
	
	def test_periodicity(self):	#transformation periodicity should be 2*pi/w_EARTH
		t1 = 1023.0
		t2 = t1 + 2*np.pi/W_EARTH
		given = np.array([5.0e3,2.30e3,5.89e3])
		v1 = fr.ecif2ecef(given,t1)
		v2 = fr.ecif2ecef(given,t2)
		self.assertTrue(np.allclose(v1,v2))
	
	@file_data('test-data/test_ecif2ecef.json')
	@unpack
	def test_cases_from_file_e2e(self,value):
		t = 12.63e3
		u = np.asarray(value)
		result = fr.ecif2ecef(u,t)
		t = t + (EPOCH - EQUINOX).total_seconds() 
		v = np.array([cos(W_EARTH*t)*value[0]+sin(W_EARTH*t)*value[1],cos(W_EARTH*t)*value[1]-sin(W_EARTH*t)*value[0],value[2]])
		self.assertTrue(np.allclose(result,v))

@ddt
class TestEcef2Ecif(unittest.TestCase):
	@file_data('test-data/test_ecif2ecef.json')
	@unpack
	def test_cases_from_file_e2e(self,value):
		t = 456.63e3
		u = np.asarray(value)
		result = fr.ecef2ecif(u,t)
		self.assertTrue(type(result),np.ndarray)
		t = t + (EPOCH - EQUINOX).total_seconds() 
		v = np.array([cos(W_EARTH*t)*value[0]-sin(W_EARTH*t)*value[1],cos(W_EARTH*t)*value[1]+sin(W_EARTH*t)*value[0],value[2]])
		self.assertTrue(np.allclose(result,v))

	def test_inverse(self):
		t = 53e2
		u = np.array([-1.0,-5.0,5e6])
		v = fr.ecef2ecif(fr.ecif2ecef(u,t),t)
		self.assertTrue(np.allclose(u,v))

@ddt
class TestEcif2Orbit(unittest.TestCase):
	@file_data('test-data/test_position_velocity.json')
	@unpack
	def test_e2o_matrix_orthogonality(self,value):
		r = np.asarray(value[0:3])
		v = np.asarray(value[3:6])
		v_x = fr.ecif2orbit(r,v,np.array([1.0,0.0,0.0]))
		v_y = fr.ecif2orbit(r,v,np.array([0.0,1.0,0.0]))
		v_z = fr.ecif2orbit(r,v,np.array([0.0,0.0,1.0]))
		self.assertEqual(type(v_z),np.ndarray)

		m = np.array([v_x,v_y,v_z])
		mT = m.transpose()
		I = np.dot(m,mT)
		self.assertTrue(np.allclose(I[0,:],[1.,0.,0.]))
		self.assertTrue(np.allclose(I[1,:],[0.,1.,0.]))
		self.assertTrue(np.allclose(I[2,:],[0.,0.,1.]))

	def test_check_pos_vel_transformation(self):
		r = np.array([1e6,-2.03,-3.0])
		v = np.array([2.0e3,2.8,-73.2])
		self.assertTrue(np.allclose(fr.ecif2orbit(r,v,r),[0.,0.,-1*np.linalg.norm(r)]))
		self.assertTrue(np.allclose(fr.ecif2orbit(r,v,np.cross(v,r)),[0.,np.linalg.norm(np.cross(v,r)),0.]))

class Test_qBI_2qBO(unittest.TestCase):	
	
	def test_south_pole_ideal(self):
		#above south pole, body frame, orbit frame and eci frame have same attitude wrt each other 
		r = np.array([0.,0.,-7.07e6])
		v = np.array([7.0e3,0.,0.])
		qBI = np.array([0.,0.,0.,1.])
		expected = np.array([0.,0.,0.,1.])
		result = fr.qBI2qBO(qBI,r,v)

		self.assertTrue(np.allclose(expected ,result))
	
	def test_data(self):
		qIB = 0.5 * np.array([0,-1,0,np.sqrt(3)])
		r = np.array([0.,0.,7.07e6])
		v = np.array([0.,7.0e3,0.])
		expectedInv = (0.5/np.sqrt(2.))*np.array([np.sqrt(3.),np.sqrt(3.),1.,1.])
		result = fr.qBI2qBO(qnv.quatInv(qIB),r,v)
		
		self.assertTrue(np.allclose(qnv.quatInv(expectedInv) ,result))

class Test_qBO_2qBI(unittest.TestCase):	

	def test_south_pole_ideal(self):
		#above south pole, body frame, orbit frame and eci frame have same attitude wrt each other 
		r = np.array([0.,0.,-7.07e6])
		v = np.array([7.0e3,0.,0.])
		qBI = np.array([1.,0.,0.,0.])
		qBO = np.array([1.,0.,0.,0.])
		result = fr.qBO2qBI(qBO,r,v)
		self.assertTrue(np.allclose(qBI ,result))	

	def test_data2(self):
		qOB = (0.5/np.sqrt(2.))*np.array([np.sqrt(3.),np.sqrt(3.),1.,1.])
		r = np.array([0.,0.,7.07e6])
		v = np.array([0.,7.0e3,0.])
		qIB = 0.5 * np.array([0,-1,0,np.sqrt(3)])
		result = fr.qBO2qBI(qnv.quatInv(qOB),r,v)		
		self.assertTrue(np.allclose(qnv.quatInv(qIB) ,result))
	
	def test_data3(self):
		r = 1/3 * np.array([7.07e6,7.07e6,7.07e6])
		v = np.array([0.,1/2**(1/2),1/2**(1/2)])
		qOB = (0.5/np.sqrt(2.))*np.array([np.sqrt(3.),np.sqrt(3.),1.,1.])
		qBI = np.array([ -0.32227667, 0.74698487, 0.43819357, 0.38227967])
		#m_OI = -1 * np.array(((2**(1/2)/9, -1/(9*2**(1/2)), -1/(9*2**(1/2))),(0, -1/(3*2**(1/2)), 1/(3*2**(1/2))),(1/3**(1/2), 1/3**(1/2), 1/3**(1/2))))
		#q_OI = qnv.quatMultiplyNorm(qnv.quatInv(qOB),qnv.rotm2quat(m_A)) = [-0.1159169  -0.88047624  0.3647052   0.27984814]

		result = fr.qBO2qBI(qnv.quatInv(qOB),r,v)	
		print(result)	
		self.assertTrue(np.allclose(qBI ,result))

	
@ddt
class Test_wBIB_2_wBOB(unittest.TestCase):

	def test_ideal_controlled(self):
		v_w_IO_o = np.array([0.,0.00106178,0.])
		qBO = np.array([1.,0.,0.,0.])
		w_BIB = -qnv.quatRotate(qBO,v_w_IO_o)
		w_BOB = fr.wBIb2wBOb(w_BIB,qBO,v_w_IO_o)
		expected = np.array([0.,0.,0.])
		self.assertTrue(np.allclose(w_BOB,expected))

	@file_data('test-data/test_stationary_body.json')
	def test_stationary_body(self,value):
		v_w_IO_o = np.array([0.,0.00106178,0.])
		qBO = np.asarray(value)
		w_BIB = np.array([0.,0.,0.])
		result = fr.wBIb2wBOb(w_BIB,qBO,v_w_IO_o)
		expected = qnv.quatRotate(qBO,v_w_IO_o)
		self.assertTrue(np.allclose(result,expected)) 

@ddt		
class Test_wBOB_2_wBIB(unittest.TestCase): #These test cases are made by reversing the input and output of Test_wBIB_2_wBOB
	def test_ideal_controlled(self): 
		v_w_IO_o = np.array([0.,0.00106178,0.])
		qBO = np.array([1.,0.,0.,0.])
		w_BOB = np.array([0.,0.,0.])
		w_BIB = fr.wBOb2wBIb(w_BOB,qBO,v_w_IO_o)
		expected = -qnv.quatRotate(qBO,v_w_IO_o)
		self.assertTrue(np.allclose(w_BIB,expected))

	@file_data('test-data/test_stationary_body.json')
	def test_stationary_body(self,value):
		v_w_IO_o = np.array([0.,0.00106178,0.])
		qBO = np.asarray(value)
		w_BOB = qnv.quatRotate(qBO,v_w_IO_o)
		result = fr.wBOb2wBIb(w_BOB,qBO,v_w_IO_o)
		expected = np.array([0.,0.,0.])
		self.assertTrue(np.allclose(result,expected)) 

if __name__=='__main__':
	unittest.main(verbosity=2)
