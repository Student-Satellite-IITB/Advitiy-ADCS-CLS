import numpy as np
import unittest
import satellite
import frames 
import default_blocks as defblock
from ddt import ddt, data, unpack, file_data
from constants_1U import v_w_IO_o, MAG_BIAS
import qnv
import sensor
import frames as fs



#Before running these test cases multiply covariance with 0 in constants_1U.py



@ddt
class TestSensor(unittest.TestCase):
	
	def test_gyroscope(self):
		qBO = np.array([.0,0.,0.,1.])
		v_w_BI_b = np.array([-3.9999, 4.8575, 0])
		v_w_BO_b = v_w_BI_b + qnv.quatRotate(qBO,v_w_IO_o)
		state = np.hstack((qBO,v_w_BO_b))
		sat = satellite.Satellite(state,12.05)
		sat.setGyroVarBias(np.array((-0.1, 0.1, 0.1)))
		self.assertTrue(np.allclose(sensor.gyroscope(sat),np.array((-4.0999,4.9575,0.1))))
	
	def test_magnetometer(self):
		qBO = np.array([.0,0.,0.,1.])
		v_w_BI_b = np.array([-3.9999, 4.8575, 0]) 
		v_w_BO_b = v_w_BI_b + qnv.quatRotate(qBO,v_w_IO_o)
		state = np.hstack((qBO,v_w_BO_b))
		sat = satellite.Satellite(state,12.05)
		v_pos_i = np.array([1e6,-2.03,-3.0])
		v_vel_i = np.array([2.0e3,2.8,-73.2])
		sat.setPos(v_pos_i)
		sat.setVel(v_vel_i)
		qOI = np.array([-0.70697614, -0.01353643,  0.70697824,  0.01353791])
		Magi = qnv.quatRotate(qnv.quatInv(qOI),np.array((2.3e-6,3.4e-6,3.1e-6)))
		sat.setMag_i(Magi) 
		self.assertTrue(np.allclose(sensor.magnetometer(sat),np.array((2.3e-6,3.4e-6,3.1e-6))))

	def test_GPS(self):
		qBO = np.array([.0,0.,0.,1.])
		v_w_BI_b = np.array([-3.9999, 4.8575, 0]) 
		v_w_BO_b = v_w_BI_b + qnv.quatRotate(qBO,v_w_IO_o)
		state = np.hstack((qBO,v_w_BO_b))
		sat = satellite.Satellite(state,3.132)
		v_pos_i = np.array([5000,5000,6000])
		v_vel_i = np.array([7e3,0,0])
		sat.setPos(v_pos_i)
		sat.setVel(v_vel_i)
		result=sensor.GPS(sat)
		self.assertTrue(np.allclose(result,np.array((5000,5000,6000,7e3,0,0,3.132))))
	
	def test_ADC_1(self): #tested corresponding to SS_GAIN = 1, SS_QUANTIZER = 3 (u=0.5) ADC_BIAS = np.array([0,0,0,0,0,0])
		result = sensor.ADC(np.array((4,5,6)))
		self.assertTrue(np.allclose(result,np.array((4,0,5,0,6,0))))

	def test_ADC_2(self): #tested corresponding to SS_GAIN = 1, SS_QUANTIZER = 3, ADC_BIAS = np.array([0,0,0,0,0,0])
		result = sensor.ADC(np.array((-4,-5,-6)))
		self.assertTrue(np.allclose(result,np.array((0,4,0,5,0,6))))

	def test_light(self): #tested corresponding to SS_THRESHOLD = 0.5
	    result=sensor.light(np.array((1,2,0.5,0.4,0.3,0)))
	    expected=np.array((1,1,0,0,0,0))
	    self.assertTrue(np.allclose(result,expected))
	
	def test_calc_SV_1(self): 
	    result=sensor.calc_SV(np.array((1,0,0,2,2,0)))
	    expected=np.array((1/(9**0.5),-2/(9**0.5),2/(9**0.5)))
	    self.assertTrue(np.allclose(result,expected))
	
	def test_calc_SV_2(self): 
	    result=sensor.calc_SV(np.array((0,2,1,0,0,3)))
	    expected=np.array((-2/(14**0.5),1/(14**0.5),-3/(14**0.5)))
	    self.assertTrue(np.allclose(result,expected))

	
	def test_sunsensor(self):
		qBO = np.array([.0,0.,0.,1.])
		v_w_BI_b = np.array([-3.9999, 4.8575, 0]) 
		v_w_BO_b = v_w_BI_b + qnv.quatRotate(qBO,v_w_IO_o)
		state = np.hstack((qBO,v_w_BO_b))
		sat = satellite.Satellite(state,12.05)
		v_pos_i = np.array([1e6,-2.03,-3.0])
		v_vel_i = np.array([2.0e3,2.8,-73.2])
		sat.setPos(v_pos_i)
		sat.setVel(v_vel_i)
		qOI = np.array([-0.70697614, -0.01353643,  0.70697824,  0.01353791])
		Suni = qnv.quatRotate(qnv.quatInv(qOI),np.array((3,4,5)))
		sat.setSun_i(Suni) 
		result = sensor.sunsensor(sat)
		expected = np.array((3/(50**(0.5)),4/(50**(0.5)),5/(50**(0.5))))
		self.assertTrue(np.allclose(result,expected))
	

if __name__=='__main__':
	unittest.main(verbosity=2)