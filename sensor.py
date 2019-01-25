from constants_1U import *
import frames as fs
import numpy as np
import qnv

def ADC(sun_vector):
    
    #works as ADC, quantizes the readings of sensor, for more you can refer to documentation of this code. 
    #input: sun vector in body frame, normal vectors of each sunsensor, sunsensor gain, quantizer
    #sunsensor gain : multiplying factor that converts dot product of sunvector with normal to voltage
    #output : 6 voltages, 1 per sunsensor
    
    u=1/(SS_QUANTIZER-1)

    m_normalVectors = np.zeros([6,3]) #matrix of normal vectors of sensors
    m_normalVectors[0,:] = v_S1
    m_normalVectors[1,:] = v_S2
    m_normalVectors[2,:] = v_S3
    m_normalVectors[3,:] = v_S4
    m_normalVectors[4,:] = v_S5
    m_normalVectors[5,:] = v_S6

    v_output = np.zeros([6]) #vector of output of each sensor

    for iter in range(0,6):
        v=np.dot(sun_vector,m_normalVectors[iter,:])
        if v<0:
            v=0  #if v < 0 that means the sunvector is making obtuse angle with the sensor which means the sensor is in dark   
        v=(u)*(round(v/u))*SS_GAIN #conversion from dot product to voltage
        v_output[iter] = v

    v_output = v_output + ADC_BIAS #add error to true quantity

    return v_output

def light(ss):
    
    #input :  voltage measured by each sunsensor, threshold voltage
    #output : light boolean per sunsensor
    #calucates the flag value-light, by giving boolean for each sun sensor 
    #according to voltage thresold 
    f=np.array([0,0,0,0,0,0])
    for i in range(6):
        if ss[i]>SS_THRESHOLD:
            f[i]=1
        else:
            f[i]=0
    return f    

def calc_SV(ss):
    #input : voltage array of all sunsensors
    #calculates back the sun vector using the sensor readings

    dark = 0 #number of sunsensors that are in dark
    sat_light = 1 #light boolean for entire satellite

    v_light = light(ss) #vector which stores whether each sunsensor is in light

    v_sun_m=np.array([0,0,0])     

    for i in range(6):
       if v_light[i]==0:
           dark=dark+1

    if dark==6:
      sat_light=0
      return np.array([0,0,0]) #if all sunsensors are dark means satellite is in dark

    if sat_light==1:     
        for n in range(3): #refer readme
              m=n*2;
              if ss[m]>=ss[m+1]:
                    v_sun_m[n]=1.0*ss[m] 
              else:
                 v_sun_m[n]=-1.0*ss[m+1]
       
    return v_sun_m/np.linalg.norm(v_sun_m)       #gives the unit sun vector to be used in quest.    


def sunsensor(sat):

    v_q_BO = sat.getQ_BO()    
    v_s_o = sat.getSun_o() #obtain sunvector in orbit frame
    v_s_b = qnv.quatRotate(v_q_BO,v_s_o) #obtain sunvector in body frame
    ss_read = ADC(v_s_b) #find reading per sunsensor
    v_sun_b_m = calc_SV(ss_read) #calculate sunvector from sunsensor readings

    return v_sun_b_m


def GPS(sat):
    #model the GPS
    #input : satellite class object
    #output : measured position, velocity, time
    #to add more outputs, add them to the hstack command and define any new constants in the constants file

    v_pos = sat.getPos() 
    v_vel = sat.getVel()
    time = sat.getTime()

    #add errors to true quantities
    v_pos_m = v_pos + GPS_POS_BIAS
    v_vel_m = v_vel + GPS_VEL_BIAS
    time_m = time + GPS_TIME_BIAS

    return np.hstack([v_pos_m,v_vel_m,time_m])

def magnetometer(sat):
    #model the magnetometer
    #input : satellite class object
    #output : measured magnetic field

    v_q_BO = sat.getQ_BO() 
    v_B_o = sat.getMag_o() #obtain magnetic field in orbit frame
    v_B_b = qnv.quatRotate(v_q_BO,v_B_o) #obtain magnetic field in body frame
    v_B_m = v_B_b + MAG_BIAS #add errors
    return v_B_m
    
def gyroscope(sat):

    v_w_BIB = sat.getW_BI_b()
    v_bias_var = sat.getGyroVarBias()
    v_w_BIB_m = v_w_BIB + v_bias_var

    return v_w_BIB_m