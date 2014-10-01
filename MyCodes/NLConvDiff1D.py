#
# Linear 1D Convection-Diffusion equation solver
#
import numpy                       
import matplotlib.pyplot as plt                 
#matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

nx = 101
dx = 2./(nx-1)
nt = 25
dt = .01
c = 1.0

cfl = dx/dt
print 'cfl = ', cfl

u0 = numpy.ones(nx)
u0[.5/dx:1/dx+1] = 2
print(u0)

plt.plot(numpy.linspace(0,2,nx), u0, color='#003366', ls='--', lw=3)
plt.ylim(0,2.5)

u = u0.copy()
for n in range(1,nt):  
    un = u.copy() 
    for i in range(1,nx): 
    
        u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1])
        
plt.plot(numpy.linspace(0,2,nx), u, color='#006633', ls='--', lw=3)
plt.ylim(0,2.5);

