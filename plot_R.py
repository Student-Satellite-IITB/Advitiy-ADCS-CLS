import matplotlib.pyplot as plt
import numpy as np

m_sgp_output_temp_i = np.genfromtxt('sgp_output_PO.csv', delimiter=",")
R=np.zeros((120000))
dot=np.zeros((120000))


for i in range(0,120000):
	R[i] = np.linalg.norm(m_sgp_output_temp_i[i,1:4])
	dot[i] = np.dot(m_sgp_output_temp_i[i,1:4]/R[i],m_sgp_output_temp_i[i,4:7]/np.linalg.norm(m_sgp_output_temp_i[i,4:7]))

plt.plot(m_sgp_output_temp_i[:,0],dot)
plt.show()