import qnv	
import unittest	
import numpy as np
from ddt import ddt,file_data,unpack,data

@ddt
class Testquatrotate(unittest.TestCase):
	@file_data("test-data/test_quatRotate.json")
	@unpack
	def test_quatrotate(self,value):
		  v = np.asarray(value[0])
		  q = np.asarray(value[1])
		  vr = np.asarray(value[2])
		  A = qnv.quatRotate(q,v)

		  self.assertTrue(np.allclose(A,vr))

@ddt
class Testquatder(unittest.TestCase):
	@file_data("test-data/test_quatder.json")
	@unpack
	def test_quatDerBI(self,value):
		  w = np.asarray(value[1])
		  q = np.asarray(value[0])
		  wdot = np.asarray(value[2])
		  A = qnv.quatDerBI(q,w)
		  self.assertTrue(np.allclose(A,wdot))

@ddt
class TestquatDerBO(unittest.TestCase):
	@file_data("test-data/test_quatDerBO.json")
	@unpack
	def test_quatDerBO(self,value):
		  w1 = np.asarray(value[1])
		  q = np.asarray(value[0])
		  qdot = np.asarray(value[2])
		  A = qnv.quatDerBO(q,w1)
		  self.assertTrue(np.allclose(A,qdot))      

@ddt
class Testrotmquat(unittest.TestCase):
	@file_data("test-data/test_rotmquat.json")
	@unpack
	def test_rotm2quat(self,value):
		   
		 
			q = np.asarray(value[3])
			m1 = np.asarray(value[0])
			m2 = np.asarray(value[1])
			m3 = np.asarray(value[2])
			A = np.vstack([m1, m2, m3])
		   
			qo = qnv.rotm2quat(A)  
			self.assertTrue(np.allclose(qo,q))       

@ddt
class TestQuat2Euler(unittest.TestCase):
	def test_zero_euler(self):
		q = np.array([0.,0.,0., 1.0])
		expected = np.array([0.,0.,0.])
		result = qnv.quat2euler(q)
		self.assertTrue(np.allclose(result,expected))
	@file_data('test-data/test_quat2euler.json')
	@unpack
	def test_quat2euler(self,value):
		q = np.asarray(value[0])
		expected = np.asarray(value[1])
		result = qnv.quat2euler(q)

		self.assertTrue(np.allclose(result,expected))

@ddt
class Testquatrotm(unittest.TestCase):
	@file_data('test-data/test_quatrotm.json')
	@unpack
	def test_quat2rotm(self,value):
			q = np.asarray(value[0])
			m1 = np.asarray(value[1])
			m2 = np.asarray(value[2])
			m3 = np.asarray(value[3])
			A = qnv.quat2rotm(q)	        
			
			self.assertTrue(np.allclose(A[0,:],m1))
			self.assertTrue(np.allclose(A[1,:],m2))
			self.assertTrue(np.allclose(A[2,:],m3))

if __name__=='__main__':
	unittest.main(verbosity=1)
