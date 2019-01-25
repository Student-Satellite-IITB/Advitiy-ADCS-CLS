import numpy as np
from constants_1U import EPOCH, EQUINOX, W_EARTH, R_EARTH
from ddt import ddt, data, unpack
from math import acos, cos, sin, pi
import unittest

'''
	This is test code for checking magfield output csv file. 
'''

@ddt
class Test_magfield(unittest.TestCase):
	
	magoutputned = np.genfromtxt('mag_output_ned.csv',delimiter=",")	
	magoutputi = np.genfromtxt('mag_output_i.csv',delimiter=",")	
	LLA = np.genfromtxt('LLA.csv',delimiter=",")
	t0 = (EPOCH - EQUINOX).total_seconds()
	T = magoutputned[:,0]
	v_result = np.zeros(4)
	l = np.linspace(0,len(T),10,int)	#Sample data points from entire file
	@data(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9])
	
	def test_magfield_predata(self,value):
		theta = (90 - self.LLA[int(value)-1, 1]) * pi/180
		if (self.LLA[int(value)-1, 2] >= 0):
			phi = (self.LLA[int(value)-1, 2]) * pi/180
		else:
			phi = (self.LLA[int(value)-1, 2] + 360) * pi/180
		m_DCMn2e = np.array([[-1*cos(theta)*cos(phi), -1*sin(phi), -1*sin(theta)*cos(phi)], [-1*cos(theta)*sin(phi), cos(phi),-1*sin(theta)*sin(phi)], [sin(theta),0.,-1*cos(theta)]])
		magoutpute = np.dot(m_DCMn2e,self.magoutputned[int(value)-1,1:4])
		
		theta = W_EARTH*(self.t0 + self.magoutputned[int(value)-1,0]) #in radian
		m_DCMe2i = np.array([[cos(theta), -1*sin(theta), 0.], [sin(theta), cos(theta),0.], [0.,0.,1.]])
		self.v_result[0] = self.magoutputned[int(value)-1,0]
		self.v_result[1:4] = np.dot(m_DCMe2i,magoutpute)
		self.magoutputi[int(value)-1,:]
		self.assertTrue(np.allclose(self.magoutputi[int(value)-1,:],self.v_result))

	magoutputitest = np.genfromtxt('mag_output_i_test.csv',delimiter=",")	
	LLAtest = np.genfromtxt('LLA_test_for_magfield.csv',delimiter=",")

	@data(0,1,2,3,4,5,6,7)	
	def test_GetLLA_postdata(self,value):
		theta = (90 - self.LLAtest[int(value), 1]) * pi/180
		if (self.LLAtest[int(value), 2] >= 0):
			phi = (self.LLAtest[int(value), 2]) * pi/180
		else:
			phi = (self.LLAtest[int(value), 2] + 360) * pi/180
		m_DCMn2e = np.array([[-1*cos(theta)*cos(phi), -1*sin(phi), -1*sin(theta)*cos(phi)], [-1*cos(theta)*sin(phi), cos(phi),-1*sin(theta)*sin(phi)], [sin(theta),0.,-1*cos(theta)]])
		magoutpute = np.dot(m_DCMn2e,self.magoutputned[int(value),1:4])
		
		theta = W_EARTH*(self.t0 + self.magoutputned[int(value),0]) #in radian
		m_DCMe2i = np.array([[cos(theta), -1*sin(theta), 0.], [sin(theta), cos(theta),0.], [0.,0.,1.]])
		self.v_result[0] = self.magoutputned[int(value),0]
		self.v_result[1:4] = np.dot(m_DCMe2i,magoutpute)
		print (self.magoutputitest[int(value),:])
		print (self.v_result)
		self.assertTrue(np.allclose(self.magoutputitest[int(value),:],self.v_result))
if __name__=='__main__':
   unittest.main(verbosity=2)