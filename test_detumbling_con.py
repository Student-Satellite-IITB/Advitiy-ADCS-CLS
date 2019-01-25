import satellite as sat
import numpy as np
import detumbling_con as detcon
import unittest	
from ddt import ddt,file_data,unpack

@ddt
class TestmagMoment(unittest.TestCase):
        @file_data("test-data/test_detumbling_con.json")
        @unpack
        
        def test_Bdot(self,value):
               
                mag1=np.asarray(value[0])
                mag2=np.asarray(value[1])
                mag_exp=np.asarray(value[2])

                iniState=np.array([1.,0.,0.,0.,0.,0.,0.])
                s = sat.Satellite(iniState,0)
                s.setMag_b_m_p(mag1)
                s.setMag_b_m_c(mag2)
                detcon.cons.k = 0.004252756437794863
                magMoment = detcon.magMoment(s)
                #print magMoment
                self.assertTrue(np.allclose(magMoment, mag_exp))
                
               
if __name__=='__main__':
        unittest.main(verbosity=2)

