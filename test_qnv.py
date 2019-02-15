import numpy as np
import qnv

A=np.array([1,2,3])
B=np.array([-1,-2,-3])
C=np.array([1,-2,-3])
D=np.array([4,-6,7])
E=np.array([0,-3,2])
F=np.array([0,0,0])

q1=np.array([4,1,2,3])
q2=np.array([-4,-1,-2,3])
q3=np.array([4,1,-2,-3])
q4=np.array([-4,-1,2,-3])
q5=np.array([4,-1,-2,-3])
q6=np.array([-4,1,2,-3])
q7=np.array([0,0,2,-3])
q8=np.array([0,0,-2,-3])
q9=np.array([0,0,0,0])

flag=1 # testing for quatInv

G=qnv.quatInv(q1)
if (G == q2).all() == 0:
	flag=0

G=qnv.quatInv(q3)
if (G == q4).all() == 0:
	flag=0

G=qnv.quatInv(q5)
if (G == q6).all() == 0:
	flag=0

G=qnv.quatInv(q7)
if (G == q8).all() == 0:
	flag=0

G=qnv.quatInv(q9)
if (G == q9).all() == 0:
	flag=0

if flag == 1:
	print ("all cases passed for quatInv")
else:
	print ("error for quatInv")
 
###############################################################################
 
flag=1 # testing for quatMultiplynorm
m1 = np.array([1,-1,-1,1])
m2 = np.array([-1,1,-1,-1])
m3 = np.array([0,1,1,-1])
m4 = np.array([0,1,-1,0])
m5 = np.array([0,0,-1,1])
m6 = np.array([-4,0,0,-0])
m7 = np.array([3,-1,-1,1])
m8 = np.array([1,1,-1,-1])
m9 = m6/np.linalg.norm(m6)
m10 = m7/np.linalg.norm(m7)
m11 = m8/np.linalg.norm(m8)

G=qnv.quatMultiplyNorm(m1,m2)
if (G == m9).all() == 0:
	flag=0

G=qnv.quatMultiplyNorm(m3,m2)
if (G == m10).all() == 0:
	flag=0

G=qnv.quatMultiplyNorm(m4,m5)
if (G == m11).all() == 0:
	flag=0

if flag == 1:
	print ("all cases passed for quatMultiplynorm")
else:
	print ("error for quatMultiplynorm")

flag=1 # testing for quatMultiplynorm

G=qnv.quatMultiplyUnnorm(m1,m2)
if (G == m6).all() == 0:
	flag=0

G=qnv.quatMultiplyUnnorm(m3,m2)
if (G == m7).all() == 0:
	flag=0

G=qnv.quatMultiplyUnnorm(m4,m5)
if (G == m8).all() == 0:
	flag=0

G=qnv.quatMultiplyUnnorm(q9,q9)
if (G == q9).all() == 0:
	flag=0

if flag == 1:
	print ("all cases passed for quatMultiplyUnnorm")
else:
	print ("error for quatMultiplyUnnorm")
