# import matplotlib.pyplot as plt
import numpy as np
from constants_1U import RESISTANCE, INDUCTANCE, PWM_AMPLITUDE, PWM_FREQUENCY, CONTROL_STEP
import math
# import time


def getEdgeCurrent(v_duty_cycle, I0):   # to return an array with current at the edges, input: duty cycle, initial current
    t_p = 1/PWM_FREQUENCY   # time period
    num_cycles = int(CONTROL_STEP/t_p)  # number of cycles=n
    dt_p = np.array(v_duty_cycle * t_p)     # time for which the voltage is high per cycle
    dt_n = np.array(t_p-dt_p)  # time for which current is low
    edgeCurrentList = np.zeros((num_cycles*2+1, 3))     # edgeCurrentList has 2n+1 rows, 3 columns for the currents
    edgeCurrentList[0, :] = I0   # setting initial value of current
    edgeCurrentList[1, :] = 1/RESISTANCE * (PWM_AMPLITUDE-(PWM_AMPLITUDE-I0*RESISTANCE)*np.exp(-RESISTANCE*dt_p/INDUCTANCE))    # Setting value of current at first edge
    for i in range(1, num_cycles):
        edgeCurrentList[2 * i, :] = edgeCurrentList[2*i-1, :]*np.exp(-RESISTANCE*dt_n/INDUCTANCE)
        edgeCurrentList[(2 * i) + 1, :] = (PWM_AMPLITUDE - (PWM_AMPLITUDE - RESISTANCE*edgeCurrentList[2 * i, :])*np.exp(-RESISTANCE*dt_p/INDUCTANCE))/RESISTANCE
    edgeCurrentList[2 * num_cycles, :] = edgeCurrentList[2 * num_cycles - 1, :] * np.exp(-RESISTANCE * dt_n / INDUCTANCE)  # Setting the current for the last edge
    return edgeCurrentList


def getAnalyticalCurrent(v_duty_cycle, edgeCurrentList, t):  #gives current at a time instant t. input: duty cycle, list of currents at edges
    t_p = 1 / PWM_FREQUENCY  # time period
    dt_p = v_duty_cycle * t_p  # time for which the voltage is high per cycle
    t_mod_tp = t % t_p  # time elapsed since last rising edge
    cur_cycle = int(t / t_p)    # gives us the index of the cycle we're in
    current_t = np.zeros(3)
    for i in range(0, 3):   # t in the comments indicates the time since last edge
        if cur_cycle >= int(CONTROL_STEP/t_p):  # if t comes after the last edge, return the current at the last edge
            current_t = edgeCurrentList[cur_cycle * 2, :]
        elif t_mod_tp < dt_p[i]:  # if the instant lies before th falling edge
            current_t[i] = (PWM_AMPLITUDE - (PWM_AMPLITUDE - RESISTANCE * edgeCurrentList[cur_cycle*2, i])*math.exp(-RESISTANCE*t_mod_tp/INDUCTANCE)) / RESISTANCE  # I=1/R(V-(V-I0R)exp(-Rt/L))
        else:
            current_t[i] = edgeCurrentList[2*cur_cycle+1, i]*math.exp(-RESISTANCE*(t_mod_tp-dt_p[i])/INDUCTANCE)    # I = I0exp(-Rt/L)
    return current_t


def getCurrentList(v_duty_cycle, t, n, edgeCurrentList): # gives current values at time instants in a given array input: duty cycle, array of time instants, number of time instants, initial current
    currentList = np.zeros((n, 3))
    for i in range(0, n):
        currentList[i, :] = getAnalyticalCurrent(v_duty_cycle, edgeCurrentList, t[i])   # current[ti]
    return currentList

'''
# this was used to compare the time taken for execution, however, it remains to be seen whether the execution in the
# actual program takes a proportional amount of time
duty_cycle = np.array([0.001, 0.000001, 0.000001])
n = 40000
t_array = np.linspace(0, 0.02, n, endpoint=False)
I0 = np.array([0, 0, 0])
# currentList = getCurrentList(duty_cycle, t_array, n, I0)
# plt.plot(t_array[:], currentList[:, 2])
# plt.show()
actstart = time.time()
current = getCurrentList(duty_cycle, t_array, n, I0)
actend = time.time()
print(actend - actstart)
print(CONTROL_STEP%(1/PWM_FREQUENCY), PWM_FREQUENCY, CONTROL_STEP, 1/PWM_FREQUENCY)
'''