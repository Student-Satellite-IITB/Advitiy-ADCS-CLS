import numpy as np
from qnv import quatRotate

def actuatorTypeA(sat):
    #input: magnetic moment required in Control Step
    #output: Torque applied in body frame about COM for duration of model step
    v_magnetic_moment_b = sat.getMagmomentRequired_b()

    v_magnetic_field_i=sat.getMag_i()
    v_magnetic_field_b=quatRotate(sat.getQ_BI(),v_magnetic_field_i) #get mag field in body frame
    
    v_torque_app_b = np.cross(v_magnetic_moment_b,v_magnetic_field_b)
    
    sat.setAppTorque_b(v_torque_app_b)

def actuatorTypeB(sat):
    #input: Torque required in body frame about COM
    #output: Torque applied in body frame about COM for duration of model step
    v_magnetic_field_i=sat.getMag_i()
    v_magnetic_field_b=quatRotate(sat.getQ_BI(),v_magnetic_field_i) #get mag field in body frame 

    v_torque_control_b = sat.getControl_b()
    #below is formulae for calculating required magnetic moment from the control torque
    v_magnetic_moment_b=(1/(np.linalg.norm(v_magnetic_field_b))**2)*np.cross(v_magnetic_field_b,v_torque_control_b)
    
    v_torque_app_b = np.cross(v_magnetic_moment_b,v_magnetic_field_b)
    
    sat.setAppTorque_b(v_torque_app_b)