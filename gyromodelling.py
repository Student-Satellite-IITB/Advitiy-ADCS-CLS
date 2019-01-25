#The following code will be used for sensor modelling of the gyroscope ITG3200
import numpy as np
import random
import satellite as sat
#The mean in x, y, z direction (ZRO) is taken to be the average of the data values taken for 10 minutes
#ZRO is assumed to be constant
#All the values that are being initialized are according to the datasheet of ITG3200
#This value of std_dev is taken from datasheet
#output = actual_value + bias + bias_change_rate(ARW) + random_errors
#ARW for now has not been taken into consideration
def gyroOoutput(sat):
    ZRO_x = -2.5899514          #(degrees/sec) zero rate output in x-direction
    ZRO_y = -1.0843420          #(degrees/sec) zero rate output in y-direction
    ZRO_z = 0.32033737          #(degrees/sec) zero rate output in z-direction
    std_dev = 0.38              #(degrees/sec)
    var=(std_dev)**2
    random_error_x = random.gauss(ZRO_x , var)  #Here we are assuming the std dev of the errors to be same for all three directions
    random_error_y = random.gauss(ZRO_y , var)
    random_error_z = random.gauss(ZRO_z , var)
    w = sat.Satellite.getW_BI_b
    output_x = (w[0] + random_error_x)            #bias
    output_y = (w[1] + random_error_y)            #bias
    output_z = (w[2] + random_error_z)            #bias
    output = np.array([output_x, output_y, output_z])
    return output

