# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 19:48:59 2014

@author: meza
"""
import numpy as np

SIZE=100
STEPS=100

#----------------------------------------------------------------------
#   Calculate internal energy
#----------------------------------------------------------------------
def energy(system, N, M):
    A1  = system[N,M]
    BCN = system[bc(N-1),M]+system[bc(N+1), M]
    BCM = system[N, bc(M-1)]+system[N, bc(M+1)]
    Etot = -A1*(BCN+BCM)
    return Etot
#----------------------------------------------------------------------#
#   Check periodic boundary conditions 
#----------------------------------------------------------------------#
def bc(i):
    if i+1 > SIZE-1:
        return 0
    if i-1 < 0:
        return SIZE-1
    else:
        return i
#----------------------------------------------------------------------
#   Build the system
#----------------------------------------------------------------------
def build_system():
    system = np.random.random_integers(0,1,(SIZE,SIZE))
    system[system==0] = -1

    return system

#----------------------------------------------------------------------
#   The Main monte carlo loop
#----------------------------------------------------------------------
def main(T):
    

    system = build_system()

    for step, x in enumerate(range(STEPS)):
        M = np.random.randint(0,SIZE)
        N = np.random.randint(0,SIZE)

        E = -2. * energy(system, N, M)

        if E <= 0.:
            system[N,M] *= -1
        elif np.exp(-1./T*E) > np.random.rand():
            system[N,M] *= -1

#----------------------------------------------------------------------#
#   Run the menu for the monte carlo simulation
#----------------------------------------------------------------------#
def run():
    print '='*70
    print '\tMonte Carlo Statistics for an ising model with'
    print '\t\tperiodic boundary conditions'
    print '='*70

    print "Choose the temperature for your run (0.1-100)"
    T = float(raw_input())
    main(T)
    