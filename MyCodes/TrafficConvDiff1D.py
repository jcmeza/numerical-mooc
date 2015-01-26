#
# Traffic Convection equation solver
#
# Author: jcmeza@ucmerced.edu
# Date: 1/26/2015
import numpy                       
import matplotlib.pyplot as plt                 
#matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16


nx = 51
L = 11. # km
dx = L/(nx-1)
dt = .001
nt = 50 # 3 minutes, with dt=.001 hour

Vmax = 136 # km/hr
rhomax = 250

#
# Set initial conditions
#
x = numpy.linspace(0,L,nx)
rho0 = numpy.ones(nx)*20.
rho0[10:20] = 50
speed0 = Vmax*(1. - (rho0/rhomax))
speed0

# Plot initial condition
plt.plot(x,rho0, color='#003366', ls='--', lw=3)
plt.ylim(0,80)

rho = rho0.copy()
for n in range(1,nt):  
    rhon = rho.copy()
    rhon[0] = 20.
    wave_speed = Vmax*(1. - 2.*(rhon/rhomax))
    for i in range(1,nx):
        rho[i] = rhon[i]-wave_speed[i]*dt/dx*(rhon[i]-rhon[i-1])
        
# Plot final solution
plt.plot(numpy.linspace(0,L,nx), rho, color='#006633', ls='--', lw=3)
plt.ylim(0,80);

