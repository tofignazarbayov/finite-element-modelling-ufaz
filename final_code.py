import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

# defining physical parameters 
L      =       # length (m)
Lambda =      # thermal conductivity of the porous domain (w/m.C)
Cp     =      # specific heat for fluid (J/kg.C)
rho    =      # density of fluid (kg/m3)
rhom   =      # density of porous domain (kg/m3)
Cpm    =      # specific heat for porous domain (J/kg.C)
u      =      # velocity (m/s)
Th     =        # hot temperature (C)
Tc     =        # cold temperature (C)
Tinit  =       # initial temperature (C)
Tmax   =        # duration of the simulation (s)
Qs     =       # heat source term (w/m3)

# defining numerical parameters 
ne     =        # number of nodes
nn     = ne + 1  # number of nodes
nu     = nn - 2  # number of unknowns
ndt    =        # number of time steps

Dx     = L / ne   # sapce step
Dt     = Tmax/ndt # time step

# defining vectors and matrix

Tu     = np.zeros(nu)          # unknown temperatures
RHS    = np.zeros(nu)          # right part of equation
A      = np.zeros((nu, nu))    # matrix
Xn     = np.linspace(0, L, nn) #vector of abscissa of nodes (used for plotting)
Temp   = np.zeros(nn)          #temperature of all nodes
Temp0  = Tinit*np.ones(nn)     #temperature of all nodes

Time   = np.zeros(ndt)         # vector for the time 
TO     = np.zeros(ndt)         # vector for the temperature at the observation point
nd     = int(ne/4) + 1 


# imposing boundary conditions 

Temp[0]      = Th
Temp[nn - 1] = Tc
Temp0[0]     = Th
Temp0[nn - 1]= Tc

# Calculating the matrix 

coeftr = rhom * Cpm *Dx/ Dt              # coef transient term
coeffD = coeftr + 2 * Lambda / Dx        # coefficients of diagonal of matrix
coeffL = -Lambda / Dx - rho * Cp * u / 2 # coefficients of upper and lower diagonals of matrix
coeffU = -Lambda / Dx + rho * Cp * u / 2 # coefficients of upper and lower diagonals of matrix

for i in range(0, nu):
    A[i, i] = coeffD

for i in range(1, nu):
    A[i, i-1] = coeffL

for i in range(0, nu-1):
    A[i, i+1] = coeffU


print(Qs*Dx)
# Calculating the RHS 

Time [0] = 0
TO   [0] = Temp0  [nd]
 


for k in range(1, ndt):
        
    for i in range(0, nu): 
        RHS [i]    = coeftr * Temp0 [i+1] + Qs * Dx
        
    
    RHS[0] = RHS[0] - coeffL * Th 
    RHS[nu - 1] = RHS [nu - 1] -coeffU * Tc 

    Tu = linalg.solve(A, RHS)
    Temp = np.concatenate((Th, Tu, Tc), axis=None)

    
    Temp0 = np.copy(Temp)
    
    Time [k] = k * Dt
    TO   [k] = Temp0  [nd]

plt.figure(1)
plt.subplot(211)
plt.plot(Xn, Temp, label = "Temperature vs x")

plt.subplot(212)
plt.plot(Time, TO, label = "Temperature vs time")

