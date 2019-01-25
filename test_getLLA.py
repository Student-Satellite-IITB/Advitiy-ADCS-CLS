import numpy as np
from constants_1U import EPOCH, EQUINOX, W_EARTH, R_EARTH
from ddt import ddt, data, unpack
from math import acos, cos, sin, pi
import unittest

'''
	This is test code for checking GetLLA output csv file. 
'''

@ddt
class Test_GetLLA(unittest.TestCase):
	
	sgpoutput = np.genfromtxt('sgp_output.csv',delimiter=",")	
	LLA = np.genfromtxt('LLA.csv',delimiter=",")	
	t0 = (EPOCH - EQUINOX).total_seconds()
	T = sgpoutput[:,0]
	LLAdata = np.zeros(4)
	
	l = np.linspace(0,len(T),10,int)	#Sample  data points from entire file
	@data(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9])
	
	def test_GetLLA_predata(self,value):
		theta = W_EARTH*(self.t0 + self.sgpoutput[int(value)-1,0]) #in radian
		m_DCM = np.array([[cos(theta), sin(theta), 0.], [-1*sin(theta), cos(theta),0.], [0.,0.,1.]])
		v_x_e = np.dot(m_DCM,self.sgpoutput[int(value)-1,1:4])
		self.LLAdata[0] = self.sgpoutput[int(value)-1,0]
		if (v_x_e[2] >= 0):
			self.LLAdata[1] = acos(((v_x_e[0]**2 + v_x_e[1]**2)/(v_x_e[0]**2 + v_x_e[1]**2 + v_x_e[2]**2))**0.5) * 180./pi
		else:
			self.LLAdata[1] = -acos(((v_x_e[0]**2 + v_x_e[1]**2)/(v_x_e[0]**2 + v_x_e[1]**2 + v_x_e[2]**2))**0.5) * 180./pi		

		if (v_x_e[1]==0):
			if (v_x_e[0]>=0):
				self.LLAdata[2] = 0.0
			else:
				self.LLAdata[2] = 180.0
		else:
			if (v_x_e[1] > 0):
				self.LLAdata[2] = acos(v_x_e[0]/((v_x_e[0]**2 + v_x_e[1]**2)**0.5))*90./(pi/2)  
			else:
				self.LLAdata[2] = -acos(v_x_e[0]/((v_x_e[0]**2 + v_x_e[1]**2)**0.5))*90./(pi/2) 

		self.LLAdata[3] = np.linalg.norm(self.sgpoutput[int(value)-1,1:4]) - R_EARTH
		v_result = self.LLA[int(value)-1,:]
		
		self.assertTrue(np.allclose(self.LLAdata,v_result))
	
	sgpoutputTest = np.genfromtxt('sgp_output_test.csv',delimiter=",")	
	t0 = (EPOCH - EQUINOX).total_seconds()
	LLATest = np.genfromtxt('LLA_test_for_getLLA.csv',delimiter=",")	

	@data(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)	
	def test_GetLLA_postdata(self,value):
		theta = W_EARTH*(self.t0 + self.sgpoutputTest[int(value),0]) #in radian
		m_DCM = np.array([[cos(theta), sin(theta), 0.], [-1*sin(theta), cos(theta),0.], [0.,0.,1.]])
		v_x_e = np.dot(m_DCM,self.sgpoutputTest[int(value),1:4])
		self.LLAdata[0] = self.sgpoutputTest[int(value),0]
		if (v_x_e[2] >= 0):
			self.LLAdata[1] = acos(((v_x_e[0]**2 + v_x_e[1]**2)/(v_x_e[0]**2 + v_x_e[1]**2 + v_x_e[2]**2))**0.5) * 180./pi
		else:
			self.LLAdata[1] = -acos(((v_x_e[0]**2 + v_x_e[1]**2)/(v_x_e[0]**2 + v_x_e[1]**2 + v_x_e[2]**2))**0.5) * 180./pi		

		if (v_x_e[1]==0):
			if (v_x_e[0]>=0):
				self.LLAdata[2] = 0.0
			else:
				self.LLAdata[2] = 180.0
		else:
			if (v_x_e[1] > 0):
				self.LLAdata[2] = acos(v_x_e[0]/((v_x_e[0]**2 + v_x_e[1]**2)**0.5))*90./(pi/2)  
			else:
				self.LLAdata[2] = -acos(v_x_e[0]/((v_x_e[0]**2 + v_x_e[1]**2)**0.5))*90./(pi/2) 

		self.LLAdata[3] = np.linalg.norm(self.sgpoutputTest[int(value),1:4]) - R_EARTH
		v_result = self.LLATest[int(value),:]
		print(self.LLAdata)
		print(v_result)
		self.assertTrue(np.allclose(self.LLAdata,v_result))

if __name__=='__main__':
   unittest.main(verbosity=2)