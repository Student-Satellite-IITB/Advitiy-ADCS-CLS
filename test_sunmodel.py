import numpy as np
from constants_1U import EPOCH, EQUINOX
from ddt import ddt, data, unpack
import unittest

'''
	This is test code for checking sun model output csv file. 
'''

@ddt
class Test_Sun_Model(unittest.TestCase):
	si_output = np.genfromtxt('si_output.csv',delimiter=",")	
	t0 = (EPOCH - EQUINOX).total_seconds()
	T = si_output[:,0]

	def calculate_sun_vector(self,t):	#calculate sun vector using code of model 
		epsilon = 23.45*np.pi/180.
		theta = (2*np.pi*(self.t0 + t)/86400.) / 365.256363
		v_sun_i = np.array([np.cos(theta), np.sin(theta)*np.cos(epsilon), np.sin(theta)*np.sin(epsilon)])
		return v_sun_i
	
	l = np.linspace(0,len(T),12,int)	#Sample  data points from entire file
	@data(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9])
	
	def test_sun_data(self,value):

		v_expected = self.calculate_sun_vector(self.T[value])
		v_result = self.si_output[value,1:4]
		
		self.assertTrue(np.allclose(v_expected,v_result))

if __name__=='__main__':
   unittest.main(verbosity=2)