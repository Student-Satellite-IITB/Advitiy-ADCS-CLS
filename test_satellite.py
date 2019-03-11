import satellite as Sat
import numpy as np
import qnv
import frames as fs
from constants_1U import v_w_IO_o

state1=np.array([1,0,0,0,10,0,0,0]) 
state2=np.array([2,0,0,0,20,0,0,0])
state3=np.array([1.1,0.1,0.1,0.1,10,0,0,0])
state4=np.array([1.2,0.2,0.2,0.2,20,0,0,0])
q=np.array([1,2,3,4])
q1 = np.array([0.5,0.5,-0.5,-0.5])
w=np.array([1,2,3])
arr1=np.array([1,2,3])
arr2=np.array([11,12,13])
arr3=np.array([1.1,2.1,3.1])
arr4=np.array([11.1,21.1,31.1])

#Case1: initializing with integer and reseting to integer

Sat1 = Sat.Satellite(state1,0)

a=Sat1.getState()

if (a == state1).all():
	print ("_init_ and getState correct")
else:
	print ("_init_ or getState incorrect")

t=Sat1.getTime()
if t==0:
	print ("_init_ and getTime correct")
else:
	print ("_init_ or getTime incorrect")

t=Sat1.setTime(1)
t=Sat1.getTime()
if t == 1:
	print ("setTime and getTime correct")
else:
	print ("setTime or getTime incorrect")

Sat1.setState(state2)
a=Sat1.getState()

if (a == state2).all():
	print ("setState correct")
else:
	print ("setState incorrect")

#Case2: initializing with integer and reseting to float
Sat2 = Sat.Satellite(state1,0)

t=Sat2.setTime(0.1)
t=Sat2.getTime()
if t == 0.1:
	print ("setTime or getTime correct")
else:
	print ("setTime or getTime incorrect")
Sat2.setState(state3)
a=Sat2.getState()

if (a == state3).all():
	print ("setState correct")
else:
	print ("setState incorrect")

#Case3: initializing with float and reseting to integer

Sat3 = Sat.Satellite(state3,0.1)

a=Sat3.getState()

if (a == state3).all():
	print ("_init_ and getState correct")
else:
	print ("_init_ or getState incorrect")

t=Sat3.getTime()
if t==0.1:
	print ("_init_ and getTime correct")
else:
	print ("_init_ or getTime incorrect")

t=Sat3.setTime(1)
t=Sat3.getTime()
if t == 1:
	print ("setTime and getTime correct")
else:
	print ("setTime or getTime incorrect")
Sat3.setState(state1)
a=Sat3.getState()

if (a == state1).all():
	print ("setState correct")
else:
	print ("setState incorrect")

#Case4: initializing with float and reseting to float
Sat4 = Sat.Satellite(state3,0.1)
t=Sat4.setTime(0.2)
t=Sat4.getTime()
if t == 0.2:
	print ("setTime correct")
else:
	print ("setTime incorrect")
Sat4.setState(state4)
a=Sat4.getState()

if (a == state4).all():
	print ("setState correct")
else:
	print ("setState incorrect")

Sat4.setPos(state1)
v_Pos_i=Sat4.getPos()

if (v_Pos_i == state1).all():
	print ("setPos or getPos correct")
else:
	print ("setPos or getPos incorrect")

Sat4.setVel(state2)
v_vel_i=Sat4.getVel()

if (v_vel_i == state2).all():
	print ("setVel or getVel correct")
else:
	print ("setVel or getVel incorrect")

Sat4.setQ_BO(q1)
Q=Sat4.getQ_BO()

if (Q==q1).all():
	print ("setQ_BO and getQ_BO correct")
else:
	print ("setQ_BO and getQ_BO incorrect")

Sat4.setW_BO_b(w)
a=Sat4.getW_BO_b()
b = Sat4.getW_BI_b()

w2 = fs.wBOb2wBIb(w,q1,v_w_IO_o)
if (a==w).all():
	print ("setW_BO_b and getW_BO_b correct")
else:
	print ("setW_BO_b and getW_BO_b incorrect")
if (b==w2).all():
	print ("getW_BI_b correct")
else:
	print ("getW_BI_b incorrect")

expected = np.array([-0.5*np.sqrt(1+0.5*np.sqrt(3.)), 
						-0.25*np.sqrt(2/(2+np.sqrt(3))), 0.25*np.sqrt(2/(2+np.sqrt(3))), 0.5*np.sqrt(1+0.5*np.sqrt(3.))])
		
Sat4.setPos(np.array([0.,0.,7.07e6]))
Sat4.setVel(np.array([0.,7.0e3,0.]))
Sat4.setQ_BO((0.5/np.sqrt(2.))*np.array([np.sqrt(3.),1.,1.,np.sqrt(3.)]))
Q=Sat4.getQ_BI()
if (np.allclose(Q, expected))==1:
	print ("getQ_BI correct")
else:
	print ("getQ_BI incorrect")

Sat4.setDisturbance_b(arr1)
a=Sat4.getDisturbance_b()

if (a==arr1).all():
	print ("setDisturbance_b and getDisturbance_b correct")
else:
	print ("setDisturbance_b or getDisturbance_b incorrect")

Sat4.setsolarDisturbance_b(arr2)
a=Sat4.getsolarDisturbance_b()

if (a==arr2).all():
	print ("setsolarDisturbance_b and getsolarDisturbance_b correct")
else:
	print ("setsolarDisturbance_b or getsolarDisturbance_b incorrect")

Sat4.setaeroDisturbance_b(arr3)
a=Sat4.getaeroDisturbance_b()

if (a==arr3).all():
	print ("setaeroDisturbance_b and getaeroDisturbance_b correct")
else:
	print ("setaeroDisturbance_b or getaeroDisturbance_b incorrect")

Sat4.setggDisturbance_b(arr4)
a=Sat4.getggDisturbance_b()

if (a==arr4).all():
	print ("setggDisturbance_b and getggDisturbance_b correct")
else:
	print ("setggDisturbance_b or getggDisturbance_b incorrect")

Sat4.setControl_b(arr1)
Sat4.setAppControl_b(arr2)

a=Sat4.getControl_b()
if (a==arr1).all():
	print ("setControl and getControl correct")
else:
	print ("setControl and getControl incorrect")

a=Sat4.getAppControl_b()
if (a==arr2).all():
	print ("setAppControl and getAppControl correct")
else:
	print ("setAppControl and getAppControl incorrect")

Sat4.setSun_i(w)
sv_i=Sat4.getSun_i()
if (sv_i==w).all():
	print ("setSun_i and getSun_i correct")
else:
	print ("setSun_i and getSun_i incorrect")

Sat4.setMag_i(arr3)
mag_i=Sat4.getMag_i()
if (mag_i==arr3).all():
	print ("setMag_i and getMag_i correct")
else:
	print ("setMag_i and getMag_i incorrect")

#Case 1 velocity and position are integer

Sat4.setPos(arr1)
v_Pos_i=Sat4.getPos()
Sat4.setVel(arr2)
v_vel_i=Sat4.getVel()

a=Sat4.getSun_o()
b=fs.ecif2orbit(v_Pos_i,v_vel_i,w)
if (a==b).all():
	print ("getSun_o correct")
else:
	print ("getSun_o incorrect")

#Case 2 velocity integer and position float

Sat4.setPos(arr3)
v_Pos_i=Sat4.getPos()
Sat4.setVel(arr2)
v_vel_i=Sat4.getVel()

a=Sat4.getSun_o()
b=fs.ecif2orbit(v_Pos_i,v_vel_i,w)
if (a==b).all():
	print ("getSun_o correct")
else:
	print ("getSun_o incorrect")

#Case 3 velocity float and position integer

Sat4.setPos(arr1)
v_Pos_i=Sat4.getPos()
Sat4.setVel(arr3)
v_vel_i=Sat4.getVel()

a=Sat4.getSun_o()
b=fs.ecif2orbit(v_Pos_i,v_vel_i,w)
if (a==b).all():
	print ("getSun_o correct")
else:
	print ("getSun_o incorrect")

#Case 4 velocity and position are float

Sat4.setPos(arr3)
v_Pos_i=Sat4.getPos()
Sat4.setVel(arr4)
v_vel_i=Sat4.getVel()

a=Sat4.getSun_o()
b=fs.ecif2orbit(v_Pos_i,v_vel_i,w)
if (a==b).all():
	print ("getSun_o correct")
else:
	print ("getSun_o incorrect")

a=Sat4.getMag_o()
b=fs.ecif2orbit(v_Pos_i,v_vel_i,mag_i)
if (a==b).all():
	print ("getMag_o correct")
else:
	print ("getMag_o incorrect")

Sat4.setSun_b_m(arr2)
a=Sat4.getSun_b_m()
if (a==arr2).all():
	print ("setSun_b_m and getSun_b_m correct")
else:
	print ("setSun_b_m and getSun_b_m incorrect")

Sat4.setMag_b_m_c(arr3)
a=Sat4.getMag_b_m_c()
if (a==arr3).all():
	print ("setMag_b_m_c and getMag_b_m_c correct")
else:
	print ("setMag_b_m_c and getMag_b_m_c incorrect")


Sat4.setMag_b_m_p(arr4)
a=Sat4.getMag_b_m_p()
if (a==arr4).all():
	print ("setMag_b_m_p and getMag_b_m_p correct")
else:
	print ("setMag_b_m_p and getMag_b_m_p incorrect")

Sat4.setQUEST(state4)
a=Sat4.getQUEST()
if (a==state4).all():
	print ("setQuest and getQuest correct")
else:
	print ("setQuest and getQuest incorrect")

Sat4.setOmega_m(state3)
a=Sat4.getOmega_m()
if (a==state3).all():
	print ("setOmega_m and getOmega_m correct")
else:
	print ("setOmega_m and getOmega_m incorrect")

Sat4.setLight(1)
if (Sat4.getLight()==1):
	print ("setLight and getLight correct")
else:
	print ("setLight and getLight incorrect")

Magmoment=np.array([-0.02,0.0,0.05])
Sat4.setMagmomentRequired_b(Magmoment)
tg=Sat4.getMagmomentRequired_b()
if (tg==Magmoment).all():
	print ("setMagmomentRequired_b and getMagmomentRequired_b correct")
else:
	print ("setMagmomentRequired_b or getMagmomentRequired_b incorrect")

Torque=np.array([-0.2,-0.03,0.5])
Sat4.setAppTorque_b(Torque)
tg=Sat4.getAppTorque_b()
if (tg==Torque).all():
	print ("setAppTorque_b and getAppTorque_b correct")
else:
	print ("setAppTorque_b and getAppTorque_b incorrect")

v_gyro_bias = np.array([-0.00001,0.0,0.3253])
Sat4.setGyroVarBias(v_gyro_bias)
bias=Sat4.getGyroVarBias()

if (bias==v_gyro_bias).all():
	print ("setGyroVarBias and getGyroVarBias correct")
else:
	print ("setGyroVarBias and getGyroVarBias incorrect")

v_gpsdata = np.array([-0.00001,0.0,0.3253,10e6,-5e6,-9,1.35])
Sat4.setgpsData(v_gpsdata)
bias=Sat4.getgpsData()
if (bias==v_gpsdata).all():
	print ("setgpsData and getgpsData correct")
else:
	print ("setgpsData or getgpsData incorrect")


v_J2data = np.array([-0.00001,0.0,0.3253,10e6,-5e6,-9,1.35])
Sat4.setJ2Data(v_J2data)
bias=Sat4.getJ2Data()
if (bias==v_J2data).all():
	print ("setJ2Data and getJ2Data correct")
else:
	print ("setJ2Data or getJ2Data incorrect")

