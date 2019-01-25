import satellite as Sat
import numpy as np
import qnv
import frames as fs

state1=np.array([1,0,0,0,10,0,0,0]) 
state2=np.array([2,0,0,0,20,0,0,0])
state3=np.array([1.1,0.1,0.1,0.1,10,0,0,0])
state4=np.array([1.2,0.2,0.2,0.2,20,0,0,0])
q=np.array([1,2,3,4])
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

Sat4.setQ(q)
Q=Sat4.getQ()

if (Q==q).all():
	print ("setQ and getQ correct")
else:
	print ("setQ and getQ incorrect")

Sat4.setW(w)
a=Sat4.getW()

if (a==w).all():
	print ("setW and getW correct")
else:
	print ("setW and getW incorrect")

Sat4.setControl_b(state1)
a=Sat4.getControl_b()
if (a==state1).all():
	print ("setControl and getControl correct")
else:
	print ("setControl and getControl incorrect")

Sat4.setSun_i(w)
sv_i=Sat4.getSun_i()
if (sv_i==w).all():
	print ("setSun_i and getSun_i correct")
else:
	print ("setSun_i and getSun_i incorrect")

Sat4.setMag_i(w)
mag_i=Sat4.getMag_i()
if (mag_i==w).all():
	print ("setMag_i and getMag_i correct")
else:
	print ("setMag_i and getMag_i incorrect")

Sat4.setSun_b_m(state2)
a=Sat4.getSun_b_m()
if (a==state2).all():
	print ("setSun_b_m and getSun_b_m correct")
else:
	print ("setSun_b_m and getSun_b_m incorrect")

Sat4.setMag_b_m(state1)
a=Sat4.getMag_b_m()
if (a==state1).all():
	print ("setMag_b_m and getMag_b_m correct")
else:
	print ("setMag_b_m and getMag_b_m incorrect")

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
if (Sat4.light==1):
	print ("setLight correct")
else:
	print ("setLight incorrect")

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


Sat4.setDisturbance_b(w)
a=Sat4.getDisturbance_b()

if (a==w).all():
	print ("setDisturbance_b and getDisturbance_b correct")
else:
	print ("setDisturbance_b or getDisturbance_b incorrect")