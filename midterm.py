import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

L = 1.0 #length
Lambda = 100
Th = 50 #hot temperature
Tc = 10 #cold temperature

ne = 5 #number of elements
nn = ne + 1 #number of nodes
nu = nn - 2 #number of unknowns
Dx = L / ne #distance between nodes

A = np.zeros((nu, nu)) #matrix
Tu = np.zeros(nu) #unknown temperatures

RHS = np.full(nu, 10) #right part of equation
Xn = np.linspace(0, L, nn) #vector of abscissa of nodes (used for plotting)
Temp = np.zeros(nn) #temperature of all nodes
Temp[0] = Th
Temp[nn - 1] = Tc

coeffD = 2 * Lambda / Dx #coefficients of diagonal of matrix
coeffLU = -Lambda / Dx #coefficients of upper and lower diagonals of matrix

for i in range(0, nu):
    A[i, i] = coeffD

for i in range(1, nu):
    A[i, i-1] = coeffLU

for i in range(0, nu-1):
    A[i, i+1] = coeffLU

RHS[0] += Lambda * Th / Dx
RHS[nu - 1] += Lambda * Tc / Dx

Tu = linalg.solve(A, RHS)

print("Matrix A:\n" + str(A))
print("RHS:\n" +str(RHS))
print("Tu:\n" + str(Tu))

Temp = np.concatenate((Th, Tu, Tc), axis=None)
plt.figure(figsize = (12, 12))
plt.plot(Xn, Temp, label = "Temperature")
plt.scatter(Xn, Temp, label = "Temperature", color = "orange")
plt.show()

