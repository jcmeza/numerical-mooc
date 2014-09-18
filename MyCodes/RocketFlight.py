import numpy as np
import matplotlib.pyplot as plt
#matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# model parameters
m_s =  50.     # weight of the rocket, kg
g = 9.8      # gravity in m s^{-2} 
rho = 1.091  # average air density, kg/m^3
r   = 0.5    # radius of rocket
A   = 2*np.pi*r**2 # cross sectional area of rocket
v_e = 325.   # exhaust speed
C_D = 0.15   # drag coefficient --- or D/L if C_L=1
m_p0 = 100   # initial weight of rocket propellant, kg

### set initial conditions ###
v0 = 0.     # start at standstill
h0 = 0.0    # initial altitude

def f(u,tn):
    """Returns the right-hand side of ODE we want to solve
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """
    
    h   = u[0]
    v   = u[1]
    m_p = np.max([m_p0 - 20.*tn,0])
    total_mass = m_s + m_p
    
    if tn <= 5.0:
        mpdot = 20.
    else:
        mpdot = 0.
        
    return np.array([v,-g + v_e/(total_mass)*mpdot - (0.5*rho*v*np.abs(v)*A*C_D)/(total_mass)])
                      
def euler_step(u, f, dt, n):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    tn = t[n]
    return u + dt * f(u,tn)


T  = 40.0                          # final time
dt = 0.1                           # time increment
N  = int(T/dt) + 1                  # number of time-steps
t  = np.linspace(0.0, T, N)      # time discretization

# initialize the array containing the solution for each time-step
u = np.empty((N, 2))
u[0] = np.array([h0, v0])# fill 1st element with initial values

# time loop - Euler method
for n in range(N-1):
    
    u[n+1] = euler_step(u[n], f, dt, n)

