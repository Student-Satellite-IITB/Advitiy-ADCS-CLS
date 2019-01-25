
import numpy as np
import math
from  constants_1U import No_Turns,v_A_Torquer,RESISTANCE,INDUCTANCE,CONTROL_STEP,h
from qnv import quatRotate

def ctrlTorqueToVoltage(sat):
    '''
        This function takes control torque from control law(this case b_dot law) and convert to magnetic moment,
        which in turn converted to current and then voltage to be applied
        Input: satellite
        Output: Voltage to be applied to Torquer to achieve the required control torque
    '''
    v_magnetic_field_i=sat.getMag_i()
    v_magnetic_field_b=quatRotate(sat.getQ(),v_magnetic_field_i) #get mag field in body frame 

    v_torque_control_b = sat.getControl_b()
    #below is formulae for calculating required magnetic moment from the control torque
    v_magnetic_moment_b=(1/(np.linalg.norm(v_magnetic_field_b))**2)*np.cross(v_magnetic_field_b,v_torque_control_b)

    Torquer_area=v_A_Torquer
    v_current=(1.0/No_Turns)*np.divide(v_magnetic_moment_b,Torquer_area)
    voltage=v_current*RESISTANCE  #simple I*R is used, since there's no other way as of now
    return voltage

def I(voltage):
    '''
        This function calculates array of current for an applied DC voltage to LR circuit
        Input: voltage to be applied as calcuated by above function(ctrlTorqueToVoltage)
        Output: array of current for a CONTROL_STEP when a voltage is directly applied to the circuit.
    '''
    N=int(CONTROL_STEP/h) #current is sampled at these many points
    t = np.linspace(0,CONTROL_STEP,N, endpoint=False)
    m_i_app_dc = np.zeros((N,4))
    for s in range(0,N):
        m_i_app_dc[s,0]=t[s]
        m_i_app_dc[s,1:4]=(voltage/RESISTANCE)*(1-np.exp(-RESISTANCE*t[s]/INDUCTANCE)) #LR circuit equation
    return  m_i_app_dc

def currentToTorque(current_list,sat):
    '''
        This function calculates the torques acting on satellite to a corresponding current in the torquer.
        Input: array of currents with first row as time and next three as currents, satellite
        Output: The torque acting on the satellite due to current in torquer.
    '''
    v_mu_app = No_Turns*np.multiply(v_A_Torquer,current_list[:,1:])
    v_magnetic_field_i=sat.getMag_i()
    v_magnetic_field_b=quatRotate(sat.getQ(),v_magnetic_field_i) #get mag field in boyd frame
    v_torque_app_b = np.cross(v_mu_app,v_magnetic_field_b)
    return v_torque_app_b
    
    

    
    
    
