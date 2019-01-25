import numpy as np
from constants_1U import m_INERTIA, M_EARTH, G, W_EARTH, AERO_DRAG, RHO, REFLECTIVITY, r_COG_2_COM_b, SOLAR_PRESSURE, Lx
import qnv as qnv
import satellite as satellite
#Refer  to EV_disturbancemodel-report_01.pdf for detailed model and equations

def ggTorqueb(s):
    '''
    this function takes input:
    Satellite Object from which it access
        position of centre of mass of satellite in eci frame
        quarternion to convert vector in eci frame to body frame
    gives output:
      torque due to gravity gradient about centre of mass in body frame
    '''
    v_q_BI = s.getQ_BI()   #unit quaternion
    v_pos_i = s.getPos()
    v_pos_b = qnv.quatRotate(v_q_BI, v_pos_i)
    pos_norm = np.linalg.norm(v_pos_b)
    
    v_t_gg_b = 3. * M_EARTH * G * (np.cross(v_pos_b, np.dot(m_INERTIA, v_pos_b))) / (pos_norm ** 5.)
    
    s.setggDisturbance_b(v_t_gg_b)

def aeroTorqueb(s):
    '''
    this function takes input:
    Satellite Object from which it access
        velocity of COM of satellite in eci frame
        quarternion to convert a vector in eci frame to body frame
    gives output:
        torque due to air drag about COM in body frame
    '''

    v_q_BI = s.getQ_BI()
    
    v_vel_sat_i = s.getVel()
    v_pos_i = s.getPos()
    v_vel_atm_i = np.cross(np.array([0.,0.,W_EARTH]),np.array([v_pos_i[0],v_pos_i[1],0.]))    #velocity of atmosphere in eci
    v_vel_i = v_vel_sat_i - v_vel_atm_i #velocity of satellite wrt atmosphere in eci 
    
    v_vel_b = qnv.quatRotate(v_q_BI, v_vel_i)
    
    area = Lx * Lx * (abs(v_vel_b[0]) + abs(v_vel_b[1]) + abs(v_vel_b[2]))  # area of satellite perpendicular to velocity
    v_t_ad_b = np.cross(r_COG_2_COM_b, v_vel_b) * RHO * AERO_DRAG  * area / 2.
    s.setaeroDisturbance_b(v_t_ad_b)

def solarTorqueb(s):
    '''
    this function takes input:
    Satellite Object from which it access
        sun vector in eci frame
        quarternion to convert a vector in eci frame to body frame
        
    gives output:
     torque due to solar drag about COM in body frame
    '''
    if s.getLight() == 0:   #in eclipse: zero solar torque
        v_t_sd_b = np.array([0.,0.,0.])
    
    else:   #light region

        v_q_BI = s.getQ_BI()
    
        v_sv_i = s.getSun_i()  # unit sun vector in inertial frame obtained from satellite object
        v_sv_b_u = qnv.quatRotate(v_q_BI, v_sv_i) / np.linalg.norm(v_sv_i)
        # area of satellite perpendicular to sun vector
        area = Lx * Lx * (abs(v_sv_b_u[0]) + abs(v_sv_b_u[1]) + abs(v_sv_b_u[2]))
        # torque due to absorption
        v_t_sd1_b = np.cross(r_COG_2_COM_b, v_sv_b_u) * SOLAR_PRESSURE * (1 - REFLECTIVITY) * area  
        # reflection torque
        v_t_sd2_b = np.cross(r_COG_2_COM_b, [abs(v_sv_b_u[0]) * v_sv_b_u[0], abs(v_sv_b_u[1]) * v_sv_b_u[1], \
                     abs(v_sv_b_u[2]) * v_sv_b_u[2]]) * 2.0* REFLECTIVITY * SOLAR_PRESSURE * Lx * Lx
        
        #total solar-pressure torque
        v_t_sd_b = v_t_sd1_b + v_t_sd2_b
        
    s.setsolarDisturbance_b(v_t_sd_b)

