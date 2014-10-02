# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 13:08:51 2014
Initial version of Burger's Equation solver
@author: meza
"""

import numpy
import sympy
import matplotlib.pyplot as plt
#%matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16
from sympy.utilities.lambdify import lambdify


from sympy import init_printing
init_printing()

x, nu, t = sympy.symbols('x nu t')
phi = sympy.exp(-(x-4*t)**2/(4*nu*(t+1))) + \
sympy.exp(-(x-4*t-2*numpy.pi)**2/(4*nu*(t+1)))
print(phi)

phiprime = phi.diff(x)
print(phiprime)

u = -2*nu*(phiprime/phi) + 4
print(u)

ufunc = lambdify((t,x,nu),u)

###variable declarations
nx = 101
nt = 100
dx = 2*numpy.pi/(nx-1)
nu = .07
dt = dx*nu

xgrid = numpy.linspace(0, 2*numpy.pi, nx)
un = numpy.empty(nx)
t = 0

u = numpy.asarray([ufunc(t, x0, nu) for x0 in xgrid])

plt.figure(figsize=(8,5), dpi=100)
plt.plot(xgrid,u, color='#003366', ls='--', lw=3)
plt.xlim([0,2*numpy.pi])
plt.ylim([0,10]);

for n in range(nt):
    un = u.copy()
    
    u[1:-1] = un[1:-1] - un[1:-1] * dt/dx * (un[1:-1] - un[:-2]) + nu*dt/dx**2*\
                    (un[2:] - 2*un[1:-1] + un[:-2])

    u[0] = un[0] - un[0] * dt/dx * (un[0] - un[-1]) + nu*dt/dx**2*\
                (un[1] - 2*un[0] + un[-1])
    u[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
                (un[0]- 2*un[-1] + un[-2])
        
u_analytical = numpy.asarray([ufunc(nt*dt, xi, nu) for xi in xgrid])

plt.figure(figsize=(8,5), dpi=100)
plt.plot(xgrid,u, color='#003366', ls='--', lw=3, label='Computational')
plt.plot(xgrid, u_analytical, label='Analytical')
plt.xlim([0,2*numpy.pi])
plt.ylim([0,10])
plt.legend();
