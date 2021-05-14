import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

L = 1 #length
Lambda = 100
Cp = 4180
rho = 1000 #density
u = 0.0001 #velocity
Th = 50 #hot temperature
Tc = 10 #cold temperature

nn = 6 #number of elements
ne = nn - 1 #number of nodes
nu = nn - 2 #number of unknowns
Dx = L / ne #distance between nodes

A = np.zeros((nu, nu)) #matrix
Tu = np.zeros(nu) #unknown temperatures

RHS = np.zeros(nu) #right part of equation

Xn = np.linspace(0, L, nn) #vector of abscissa of nodes (used for plotting)
Temp = np.zeros(nn) #temperature of all nodes
Temp[0] = Th
Temp[nn - 1] = Tc

coeffD = 2 * Lambda / Dx #coefficients of diagonal of matrix
coeffL = -Lambda / Dx - rho * Cp * u / 2 #coefficients of upper and lower diagonals of matrix
coeffU = -Lambda / Dx + rho * Cp * u / 2 #coefficients of upper and lower diagonals of matrix

for i in range(0, nu):
    A[i, i] = coeffD

for i in range(1, nu):
    A[i, i-1] = coeffL

for i in range(0, nu-1):
    A[i, i+1] = coeffU

RHS[0] = -coeffL * Th
RHS[nu - 1] = -coeffU * Tc

Tu = linalg.solve(A, RHS)

print("Matrix A:\n" + str(A))
print("RHS:\n" +str(RHS))
print("Tu:\n" + str(Tu))

Temp = np.concatenate((Th, Tu, Tc), axis=None)
plt.plot(Xn, Temp, label = "Temperature")
plt.scatter(Xn, Temp, label = "Temperature", color = "orange")
plt.show()