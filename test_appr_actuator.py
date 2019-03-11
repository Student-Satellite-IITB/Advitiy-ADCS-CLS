import satellite as sat
import numpy as np
import appr_actuator as act
from qnv import quatInv, quatRotate
from frames import qBI2qBO
import unittest 
from ddt import ddt,file_data,unpack

@ddt
class TestapprActuator(unittest.TestCase):
        @file_data("test-data/test_appr_actuatorTypeA.json")
        @unpack
        
        def test_appr_actTypeA(self,value):
               
                mag_moment_b=np.asarray(value[0])
                magfield_b=np.asarray(value[1])
                app_torque_exp=np.asarray(value[2])

                vel=np.array([5,6,7])
                pos=np.array([3,4,5])

                q_BI=np.array([0.5,-0.5,0.5,-0.5])
                q_IB=quatInv(q_BI)
                q_BO=qBI2qBO(q_BI,pos,vel)
    
                magfield_i=quatRotate(q_IB, magfield_b)
                
                w_BO_b=np.array([8,6,7])
                state=np.hstack((q_BO,w_BO_b))
                t0=0
                s=sat.Satellite(state,t0)
                s.setVel(vel)
                s.setPos(pos)

                s.setMagmomentRequired_b(mag_moment_b)
                s.setMag_i(magfield_i)
                
                app_torque=act.actuatorTypeA(s)
                
                self.assertTrue(np.allclose(app_torque, app_torque_exp))

        @file_data("test-data/test_appr_actuatorTypeB.json")
        @unpack
        
        def test_appr_actTypeB(self,value):
               
                torque_req=np.asarray(value[0])
                magfield_b=np.asarray(value[1])
                app_torque_exp=np.asarray(value[2])

                vel=np.array([5,6,7])
                pos=np.array([3,4,5])

                q_BI=np.array([0.5,-0.5,0.5,-0.5])
                q_IB=quatInv(q_BI)
                q_BO=qBI2qBO(q_BI,pos,vel)

                magfield_i=quatRotate(q_IB, magfield_b)

                w_BO_b=np.array([8,6,7])
                state=np.hstack((q_BO,w_BO_b))
                t0=0
                s=sat.Satellite(state,t0)
                s.setVel(vel)
                s.setPos(pos)
                s.setMag_i(magfield_i)
                s.setControl_b(torque_req)
                
                app_torque=act.actuatorTypeB(s)
                
                self.assertTrue(np.allclose(app_torque, app_torque_exp))
                
               
if __name__=='__main__':
        unittest.main(verbosity=2)


