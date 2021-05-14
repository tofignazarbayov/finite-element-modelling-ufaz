import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

L = 1.0 #length (m)
Lambda = 100   # thermal conductivity (w/m.C)
Cp     = 4180  # specific heat for fluid (J/kg.C)
rho    = 1000  # density of fluid  (kg/m3)
rhom   = 2000  # density of porous domain  (kg/m3)
Cpm    = 3000  # specific heat for fluid (J/kg.C)
u      = 0.01  # velocity (m/s)
Th     = 50    # hot temperature (C)
Tc     = 10    # cold temperature (C)
Tinit  = 10    # initial temperature (C)
Tmax   = 3     # duration of the simulation (s)  

nn  = 5      # number of elements
ne  = nn - 1 # number of nodes
nu  = nn - 2 # number of unknowns
ndt = 3      # number of time steps

Dx = L / ne   # sapce step
Dt = Tmax/ndt # time step

coefftr = rhom * Cpm * Dx / Dt

A = np.zeros((nu, nu)) #matrix

Tu = np.zeros(nu) #unknown temperatures
#RHS = coefftr * Tinit * np.ones(nu) #right part of equation
RHS = np.full(nu, coefftr * Tinit)

Xn = np.linspace(0, L, nn) #vector of abscissa of nodes (used for plotting)

Temp = np.zeros(nn) #temperature of all nodes
Temp[0] = Th
Temp[nn - 1] = Tc

coeffD = coefftr + 2 * Lambda / Dx #coefficients of diagonal of matrix
coeffL = -Lambda / Dx - rho * Cp * u / 2 #coefficients of upper and lower diagonals of matrix
coeffU = -Lambda / Dx + rho * Cp * u / 2  #coefficients of upper and lower diagonals of matrix

for i in range(0, nu):
    A[i, i] = coeffD

for i in range(1, nu):
    A[i, i-1] = coeffL

for i in range(0, nu-1):
    A[i, i+1] = coeffU

RHS[0] = -coeffL * Th #+ rhom * Cpm * Tinit * Dx / Dt
RHS[nu - 1] = -coeffU * Tc #+ rhom * Cpm * Tinit * Dx / Dt

Tu = linalg.solve(A, RHS)

Temp = np.concatenate((Th, Tu, Tc), axis=None)
plt.plot(Xn, Temp, label = "Temperature")
plt.scatter(Xn, Temp, label = "Temperature")
plt.show()