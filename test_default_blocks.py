import numpy as np
import qnv
import unittest
import satellite
import frames 
import default_blocks as defblock
from ddt import ddt, data, unpack, file_data

@ddt
class TestDefaultBlocks(unittest.TestCase):
	
	@file_data("test-data/test_defaultsunsensor.json")
	@unpack
	def test_sunsensor(self,value): 
		qBI = np.asarray(value[0])
		v_sv_i=np.asarray(value[1])

		v_pos_i = np.array([1e6,-2.03,-3.0])
		v_vel_i = np.array([2.0e3,2.8,-73.2])
		qBO = frames.qBI2qBO(qBI,v_pos_i,v_vel_i)
		state = np.hstack((qBO,np.zeros(3)))
		mySat = satellite.Satellite(state,12.0)
		mySat.setSun_i(v_sv_i)
		mySat.setPos(v_pos_i)
		mySat.setVel(v_vel_i)
		qOI = qnv.quatMultiplyNorm(qnv.quatInv(qBO),qBI)
		result = defblock.sunsensor(mySat)
		v_expected = value[2];

		self.assertTrue(np.allclose(result,v_expected))
	
	@file_data("test-data/test_defaultsunsensor.json")
	@unpack
	def test_magmeter(self,value): 
		qBI = np.asarray(value[0])
		v_B_i=np.asarray(value[1])
		v_pos_i = np.array([1e6,-2.03,-3.0])
		v_vel_i = np.array([2.0e3,2.8,-73.2])
		qBO = frames.qBI2qBO(qBI,v_pos_i,v_vel_i)
		state = np.hstack((qBO,np.zeros(3)))
		mySat = satellite.Satellite(state,12.0)
		mySat.setMag_i(v_B_i)
		mySat.setPos(v_pos_i)
		mySat.setVel(v_vel_i)
		qOI = qnv.quatMultiplyNorm(qnv.quatInv(qBO),qBI)
		result = defblock.magnetometer(mySat)
		v_expected = value[2];
		self.assertTrue(np.allclose(result,v_expected))
	
	def test_gyroscope(self):
		qBO = np.array([.0,0.,0.,1.])
		state = np.hstack((qBO,np.array([-0.1,0.39159385,-0.7])))
		sat = satellite.Satellite(state,12.05)
		v_w_BI_b = sat.getW_BI_b()
		self.assertTrue(np.allclose(defblock.gyroscope(sat),v_w_BI_b))

	def test_GPS(self):
		qBO = np.array([.0,0.,0.,1.])
		state = np.hstack((qBO,np.array([-0.1,0.39159385,-0.7])))
		sat = satellite.Satellite(state,12.05)
		sat.setPos(np.array([7070e3,0.,0.])) 
		sat.setVel(np.array([5.60,-5.0,0.0]))
		v_expected = np.array([7070e3,0.,0.,5.60,-5.0,0.0,12.05])
		self.assertTrue(np.allclose(defblock.gps(sat),v_expected))


	def test_J2_propagator(self):
		qBO = np.array([.0,0.,0.,1.])
		state = np.hstack((qBO,np.array([-0.1,0.39159385,-0.7])))
		sat = satellite.Satellite(state,12.05)
		sat.setPos(np.array([7070e3,0.,0.])) 
		sat.setVel(np.array([5.60,-5.0,0.0]))
		v_expected = np.array([7070e3,0.,0.,5.60,-5.0,0.0,12.05])
		self.assertTrue(np.allclose(defblock.J2_propagator(sat),v_expected))

	def test_controller(self):
		qBO = np.array([.0,0.,0.,1.])
		state = np.hstack((qBO,np.array([0.10050588,-0.05026119,-0.3014887])))
		mySat = satellite.Satellite(state,128.05)
		self.assertTrue(np.allclose(defblock.controller(mySat),np.zeros(3)))

	def test_disturbance(self):
		qBO = np.array([.0,0.,0.,1.])
		state = np.hstack((qBO,np.array([0.10050588,-0.05026119,-0.3014887])))
		mySat = satellite.Satellite(state,128.05)
		self.assertTrue(np.allclose(defblock.disturbance(mySat),np.zeros(3)))
	
	def test_estimator(self):
		qBO = np.array([.414,0.5,-0.5,0.])
		state = np.hstack((qBO,np.array([0.10050588,-0.05026119,-0.3014887])))
		mySat = satellite.Satellite(state,128.05)
		self.assertTrue(np.allclose(defblock.estimator(mySat),qBO))
	
if __name__=='__main__':
	unittest.main(verbosity=2)