import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# model parameters
m_s = 50.     # weight of the rocket, kg
g = 9.81     # gravity in m s^{-2} 
rho = 1.091  # average air density, kg/m^3
r = 0.5    # radius of rocket
A = np.pi*r**2 # cross sectional area of rocket
v_e = 325.   # exhaust speed
C_D = 0.15   # drag coefficient
m_p0 = 100.   # initial weight of rocket propellant, kg

### set initial conditions ###
v0 = 0.0     # start at standstill
h0 = 0.0    # initial altitude

def f(u,n):
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
    
    h = u[0] # not used here, but included for clarity
    v = u[1]
    
    mdot = mpdot[n]
        
    m_p = m_p0 - mdot*t[n]

    if t[n] >= 5.0:
        m_p = 0.0 # all propellant used up after 5 seconds
        
    total_mass = m_s + m_p
    
    vrhs = -g + (v_e*mdot)/total_mass -(0.5*rho*v*np.abs(v)*A*C_D)/total_mass
    
    return np.array([v, vrhs])
                      
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
    return u + dt*f(u,n)


T = 40.0                          # final time
dt = 0.1
                           # time increment
N = int(T/dt) + 1                  # number of time-steps
t = np.linspace(0.0, T, N)      # time discretization

#
mpdot = np.zeros(t.shape)
mpdot[np.where(t<=5.)] = 20.

# initialize the array containing the solution for each time-step
u = np.empty((N, 2))
u[0] = np.array([h0, v0])# fill 1st element with initial values

# time loop - Euler method
for n in range(N-1):
    u[n+1] = euler_step(u[n], f, dt, n)
    print t[n+1], u[n+1]
    
# Answers to homework questions
    
height = u[:,0]
velocity = u[:,1]
vmax = np.max(velocity)
t_at_vmax = t[np.argmax(velocity)]
h_at_vmax = height[np.argmax(velocity)]
print 'Vmax = ', vmax, 'at time = ', t_at_vmax, 'with height = ', h_at_vmax

hmax = np.max(height)
t_at_hmax = t[np.argmax(height)]
print 'Hmax = ', hmax, 'at time = ', t_at_hmax

# plot height
plt.plot(t,height)

# plot velocity
plt.plot(t,velocity)
