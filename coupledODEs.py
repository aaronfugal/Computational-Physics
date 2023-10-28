import numpy as np
from matplotlib import pyplot as plt

k = 0.00033
g = 9.8

Dt = 0.01 #time step
vy0 = 0 #initial vy
vx0 = 45 # initial vx
t_start = 0
t_end = 60

n_steps = int(round((t_end-t_start)/Dt))

vx_arr = np.zeros(n_steps +1)
vy_arr = np.zeros(n_steps +1)
t_arr = np.zeros(n_steps +1)
t_arr[0] = t_start
vx_arr[0] = vx0
vy_arr[0]= vy0

for i in range (1, n_steps +1):
    vy = vy_arr[i-1]
    vx = vx_arr[i-1]
    t = t_arr[i-1]
    Dvydt = -k*vy*np.sqrt(vx**2+vy**2) - g 
    Dvxdt= -k*vx*np.sqrt(vx**2+vy**2)
    vy_arr[i] = vy + Dt*Dvydt
    vx_arr[i] = vx + Dt*Dvxdt
    t_arr[i] = t + Dt
    
    

    


fig = plt.figure()                                  # create figure
plt.plot(t_arr, vx_arr, linewidth = 1, label = 'vx')    # plot Y to t 
plt.plot(t_arr, vy_arr, linewidth = 5, label = 'vy')    # plot P to t
plt.title('Title', fontsize = 12)    # add some title to your plot
plt.xlabel('t (in seconds)', fontsize = 12)
plt.ylabel('vx(t), vy(t)', fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)                        # show grid
plt.axis([t_start, t_end, 0, 50])     # show axes measures
plt.legend()
plt.show()

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pylab

########################
### global constants ###
########################

g     = 9.8               # m/s**2
v0x   = 45.0              # m/sec
v0y   = 0.0               # m/sec
k     = 0.0033            # m**-1
x0    = 0.0               # m
y0    = 80.0              # m
tmax  = np.sqrt(2.0*y0/g) # sec; time to fall without air resistance
dt  = 0.01              # sec
n  = int(tmax/dt)

#################
def dvxdt(vx,vy):
#################

    return -1.0*k*vx*np.sqrt(vx*vx+vy*vy)
    
#################
def dvydt(vx,vy):
#################

    return -1.0*k*vy*np.sqrt(vx*vx+vy*vy) - g
    
###################
def dxdt(vx):
###################

    return vx

###################
def dydt(vy):
###################

    return vy

###################
### main script ###
###################

t = np.arange(0.0,tmax,dt)   # time parameter

vx    = np.zeros(len(t))
vy    = np.zeros(len(t))
x     = np.zeros(len(t))
y     = np.zeros(len(t))
x[0]  = x0
y[0]  = y0
vx[0] = v0x
vy[0] = v0y

### first, calculate the velocities without air resistance ###

vxvac = np.ones(len(t))*v0x
vyvac = v0y - g*t

### now, evolve vx(t) and vy(t) using Euler method ###

for i in range(1,n+1):

    vx[i] = vx[i-1] + dt/2*(dvxdt(vx[i-1],vy[i-1])+dvxdt(vx[i],vy[i]))
    vy[i] = vy[i-1] + dt/2*(dvydt(vx[i-1],vy[i-1])+dvydt(vx[i],vy[i]))
    x[i]  = x[i-1] + dt/2*(dxdt(vx[i-1])+dxdt(vx[i]))
    y[i]  = y[i-1] + dt/2*(dydt(vy[i-1])+dydt(vy[i]))
    
xvac = vx*t
yvac = y0 + vy*t

### plot results, compare with no air resistance ###


fig = plt.figure(figsize=(8,8))

fig.add_subplot(211)
plt.plot(t,vx,   lw=2,c='b',linestyle='-',label='$v_x$')
plt.plot(t,vxvac,lw=2,c='b',linestyle=':',label='$v_x$ (vacuum)')
plt.plot(t,vy,   lw=2,c='r',linestyle='-',label='$v_y$')
plt.plot(t,vyvac,lw=2,c='r',linestyle=':',label='$v_y$ (vacuum)')
plt.xlabel('time (seconds)')
plt.ylabel('speed (m/s)')
plt.grid(True)
plt.legend(loc="lower left")
fig.add_subplot(212)
plt.plot(x,y, lw=2, c='g', label = 'air resistance')
plt.plot(xvac,yvac, lw=2, c='m', label = 'vacuum')
plt.grid(True)
plt.legend(loc = 'lower left')
plt.xlabel('x (m)')
plt.ylabel('y (m)')

plt.show()
