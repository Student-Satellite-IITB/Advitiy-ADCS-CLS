import numpy as np
from constants_1U import R_EARTH, AU, R_SUN
'''
    This code generates light model output file for given orbit
    Input: orbit data file, sun vector data file
    output: N*2 array data file 
    flag = 0 eclipse
    flag = 1 light
    flag = 0.5 penumbra
    For details of model refer: 
'''


#Read SGP and sun-model data
m_sgp_output = np.genfromtxt('sgp_output.csv', delimiter=",")
m_si_output = np.genfromtxt('si_output.csv', delimiter=",")
T = m_sgp_output[:,0] #storing first element as time

N = len(T)
m_light_output = np.zeros((N,2))

r_umbra = AU*R_EARTH / (R_SUN - R_EARTH) #distance from vertex of the umbra cone to center of the earth
r_penumbra = AU*R_EARTH / (R_SUN + R_EARTH) #distance from vertex of the penumbra cone to center of the earth
alpha = np.arcsin(R_EARTH / r_umbra) #Half the aperture of the cone made by umbra (in radians)
beta = np.arcsin(R_EARTH / r_penumbra) #Half the aperture of the cone made by penumbra (in radians)

for i in range(N):
    v_pos_i = m_sgp_output[i,1:4].copy()  #position of satellite in ECI
    v_sun_i = m_si_output[i,1:4].copy()   #sun vector in ECI
    #angle between sun-vector and satellite position vector in ECI frame
    theta = np.arccos(np.dot(v_pos_i, v_sun_i) /np.linalg.norm(v_pos_i) )
    #angle between the sunvector and vector from the vertex of the umbra cone to satellite
    theta_u = np.arccos((np.dot((v_pos_i + r_umbra*v_sun_i), v_sun_i)) / np.linalg.norm(v_pos_i + r_umbra*v_sun_i))
    #angle between the negative sunvector and vector from the vertex of the penumbra cone to satellite
    theta_p =np.arccos((np.dot((v_pos_i - r_penumbra*v_sun_i), -v_sun_i)) / np.linalg.norm(v_pos_i - r_penumbra*v_sun_i))
    
    #Boolean to store whether satellite is in light or dark. 1 implies satellite is in light.
    if (theta >= np.pi/2 + alpha) & (theta_u <= alpha):
        flag = 0
    elif (theta >= np.pi/2 - beta)  & (theta_p <= beta):
        flag = 0.5
    else:
    	flag = 1

    m_light_output[i,0] = T[i]
    m_light_output[i,1] = flag
    
np.savetxt("light_output.csv", m_light_output, delimiter=',')
